#!/usr/bin/env python
# coding: utf-8

# In this notebook, starting from diabetic parameters (from Dalla Man paper), we infer a parameters subsets (conresponding to each compartment), in order to fit data simulated from healthy subject parameters also taken from Dalla Man (instead of Pattou's data).

# computer algebra library
# import symengine
import sympy

from sympy import Eq, Derivative, Piecewise, DiracDelta, ITE, O

# numeric integration library
import scipy
from scipy.integrate import solve_ivp

# plotting library
from matplotlib.pyplot import plot, show, legend, xlabel, ylabel, title, xlim, ylim, grid, figure

# library for scientific computing, used here for plotting
import numpy

from IPython.display import Image

from customDisplays import *

from customSubstitutions import *

# serialization
import pickle
import json
import dill

import pandas as pd

from lmfit import minimize, Parameters, Model

from model import *

from parametersValues import parameter_values

# regular expressions
import re

# Loading the steady state symbolic basal values
lmodel_ss_sbv = pickle.load(open('model_ss_sbv.p', 'rb'))


# Initial conditions (*We should check the value for D*)
initial_conditions = {
    Gp: Gpb,    # eq (1)
    Gt: Gtb,    # eq (1)
    G: Gb,      # eq (1)
    EGP: EGPb,  # eq (10)
    Il: Ilb,    # eq (3)
    Ip: Ipb,    # eq (3)
    I: Ib,      # eq (3)
    I1: Ib,     # eq (11)
    Id: Ib,     # eq (11)
    HE: HEb,    # eq (4)
    Qsto: 0,    # eq (13)
    Qsto1: 0,   # eq (13)
    Qsto2: 0,   # eq (13)
    Qgut: 0,    # eq (13)
    Ra: 0,      # eq (13)
    X: 0,       # eq (18)
    Y: 0,       # eq (26)
    Ipo: Ipob,  # eq (24)
    D: 90000
}


# ## Loading serialized python objects

# symbolic equations
toexp = json.load(open('tostr.json', 'r'))
toexp = eval(toexp)

initial_odes = {}

for k, v in toexp.items():
    initial_odes[eval(k)] = eval(v)

# ## Preparing model for parameter estimation
# To sum up:
# - `lmodel_ss_sbv` : symbolic basal steady state values (maps basal symbols to
#    symbolic expressions)
# - `initial_conditions` : the initial conditions (maps variables to 0 or basal symbol)
# - `sympy_odes` : the symbolic differential equations
# - `gly_dataset` : the data set to fit

# We restore D as being a constant rather than a function symbol and substitute the Diract function with a sigmoid. Also, we substitute h for Gb in the equation for Y:

Dconst = SSymbol('Dconst')
initial_odes.update({Qsto1: initial_odes[Qsto1].subs(D, Dconst).subs(
    DiracDelta(t), (0.5 - sympy.tanh(10*(abs(t)-1))/2))})
initial_odes.update({Qsto2: initial_odes[Qsto2].subs(D, Dconst)})
initial_odes.update({Qgut: initial_odes[Qgut].subs(D, Dconst)})
initial_odes.update({Y: initial_odes[Y].subs(h, Gb)})

# `HEb` is already present in the parameter values, we thus remove it from the basal symbols:
del lmodel_ss_sbv[HEb]

healthy_values = parameter_values['normal'].copy()
diabetic_values = parameter_values['diabetic'].copy()
diabetic_values.update({Dconst: initial_conditions[D]})
healthy_values.update({Dconst: initial_conditions[D]})
parameter_names = tuple(n for n in diabetic_values.keys())
parameter_names_str = tuple(str(n) for n in parameter_names)
# diabetic_values_str maps string parameter names to symbolic expressions:
diabetic_values_str = {pns: diabetic_values[eval(pns)] for pns in parameter_names_str}
variables = tuple([v for v in initial_odes.keys()])
basal_symbols = tuple(lmodel_ss_sbv.keys())
basal_symbols_str = tuple(str(bs) for bs in basal_symbols)

# We restrict the initial conditions to the variables of interest:
initial_conditions = {v: initial_conditions[v] for v in variables}

# We first define a function that simplifies as much as possible all of the symbolic expressions:
def partial_evaluation(params, param_to_estimate):
    '''
    This function performs the partial evaluation of all the parameter_values but the parameter to
    estimate to the odes, the basal symbolic expressions and the initial conditions.
    Inputs: - parameter_values: parameter values ie. diabetic or normal
    - param_to_estimate: the parameters that should remain symbolic
    Outputs: - a dictionnary with parameter_values minus param_to_estimate, the basal symbolic
    expressions, the initial conditions and the odes
    '''
    parameter_values = params.copy()
    odes = initial_odes.copy()
    symb_basal_expr = lmodel_ss_sbv.copy()
    init_cond = initial_conditions.copy()
    # we keep in parameters only the parameter values that we don't want to estimate:
    for p in param_to_estimate:
        del parameter_values[p]
    # we substitute these parameters for their value in symbolic basal expressions:
    for bs in basal_symbols:
        symb_basal_expr[bs] = symb_basal_expr[bs].subs(parameter_values) \
                              if (type(symb_basal_expr[bs]) != int) \
                              else symb_basal_expr[bs]
    # ... in odes:
    for v in variables:
        odes[v] = odes[v].subs(symb_basal_expr).subs(parameter_values) \
                  if (type(odes[v]) != int) \
                  else odes[v]
    # ... and in initial conditions:
    for v in variables:
        init_cond[v] = init_cond[v].subs(symb_basal_expr) \
                       if (type(init_cond[v]) != int) \
                       else init_cond[v]
    return {
             "parameter_values" : parameter_values,
             "symb_basal_expr": symb_basal_expr,
             "init_cond": init_cond,
             "odes": odes
           }

# ## Generate the data to fit from healthy subject
def generate_simulated_data(params, nb_points=30):
    '''
    This function generates data to be fitted for a given set of parameter values.
    Inputs: - params: the parameter values of the model for which the data should be
    simulated.
    - nb_points: the number of points that should be generated.
    Outputs: - two arrays of data (one for plasma glucose (Gp) and one for plasma insulin (Ip)).
    '''
    parameter_values = params.copy()
    time_interv = scipy.linspace(0, 300, nb_points)
    simplified_model = partial_evaluation(parameter_values, {})
    (odes,init_cond) = (simplified_model['odes'], simplified_model['init_cond'])
    lambdified_odes = [
                        sympy.lambdify((t,) + variables, odes[v].args[1], modules='sympy')
                        for v in variables
                      ]
    def numsys(t, y, *kargs, **kwargs): return [eqn(t, *y, **kwargs) for eqn in lambdified_odes]
    sol_healthy = solve_ivp(
                            numsys,
                            [0, 300],
                            tuple(init_cond[v] for v in variables),
                            method="LSODA",
                            dense_output=True
                           )
    data = [
            sol_healthy.sol(time_interv)[0],
            sol_healthy.sol(time_interv)[variables.index(Ip)]
           ]
    return data

# plot(time_interv,data[0],'b.')
# plot(time_interv,data[1],'r.')
#legend(["glucose data", "insulin data"])


# ## Lambdifying the ODEs
def lambdify_odes(param_to_estimate, odes):
    '''
    Lambdifies the odes which are abstracted away from the variables and the parameter to estimate
    '''
    return [sympy.lambdify((t,) + variables + param_to_estimate,
                           odes[v].args[1], modules='sympy')
            for v in variables
            ]


def numsys(t, param_to_estimate, y, odes, *kargs, **kwargs):
    return [eqn(t, *y, **kwargs) for eqn in lambdify_odes(param_to_estimate, odes)]

def lambdify_init_cond(param_to_estimate, init_cond):
    '''
    Lambdifies the initial conditions w.r.t. the parameters to estimate.
    Input: - parameters to estimate
    Output: - a dictionary that associates to each variable v a lambda term mapping the parameters
              to estimate to the steady state value of v
    '''
    return { v : sympy.lambdify(param_to_estimate, init_cond[v], modules=sympy) for v in variables }

def my_function(t, param_to_estimate, lambdified_init_cond, odes, **param_to_estimate_vals_str):
    '''
    Function used by the fitting algorithms to compute the my_function error.
    '''
    param_values = tuple(param_to_estimate_vals_str[str(p)] for p in param_to_estimate)
    numerical_init_cond = [ lambdified_init_cond[v](*param_values) for v in variables ]
    sol_num = solve_ivp(lambda t, y: \
                        numsys(t, param_to_estimate, y, odes, **param_to_estimate_vals_str),
                        (0, 300),
                        numerical_init_cond,
                        method="LSODA",
                        dense_output=True)
    return numpy.array([sol_num.sol(t)[0], sol_num.sol(t)[variables.index(Ip)]])

# ## Parameter fitting for the intestinal absorption compartment
def generate_intestinal_model():
    '''
    Generates a model for the intestinal compartment.
    '''
    param_to_estimate = (kmax, kmin, kab, kgri, f, b, c, Gb)
    param_to_estimate_str = [str(p) for p in param_to_estimate]
    symplified_model = partial_evaluation(diabetic_values, param_to_estimate)
    params = symplified_model["parameter_values"]
    symb_basal_expr = symplified_model["symb_basal_expr"]
    init_cond = symplified_model["init_cond"]
    odes = symplified_model["odes"]
    lambdified_init_cond = lambdify_init_cond(param_to_estimate, init_cond)
    return Model(lambda t, p_kmax, p_kmin, p_kab, p_kgri, p_f, p_b, p_c, p_Gb :
                    my_function(t,
                                param_to_estimate,
                                lambdified_init_cond,
                                odes,
                                **{"kmax":p_kmax, "kmin":p_kmin, "kab":p_kab, "kgri":p_kgri,
                                   "f":p_f, "b":p_b, "c":p_c, "Gb":p_Gb }))

# ## Parameter fitting for the glucose kinetics compartment
def generate_glucose_kinetics_model():
    '''
    Generates a model for the intestinal compartment.
    '''
    param_to_estimate = (VG, k1, k2, Gb)
    param_to_estimate_str = [str(p) for p in param_to_estimate]
    symplified_model = partial_evaluation(diabetic_values, param_to_estimate)
    params = symplified_model["parameter_values"]
    symb_basal_expr = symplified_model["symb_basal_expr"]
    init_cond = symplified_model["init_cond"]
    odes = symplified_model["odes"]
    lambdified_init_cond = lambdify_init_cond(param_to_estimate, init_cond)
    return Model(lambda t, p_VG, p_k1, p_k2, p_Gb :
                    my_function(t,
                                param_to_estimate,
                                lambdified_init_cond,
                                odes,
                                **{ "VG":p_VG, "k1":p_k1, "k2":p_k2, "Gb": p_Gb }))

# ## Parameter fitting for the insulin kinetics compartment
def generate_insulin_kinetics_model():
    '''
    Generates a model for the intestinal compartment.
    '''
    param_to_estimate = (VI, m1, m2, m4, m5, m6, HEb, Gb)
    param_to_estimate_str = [str(p) for p in param_to_estimate]
    symplified_model = partial_evaluation(diabetic_values, param_to_estimate)
    params = symplified_model["parameter_values"]
    symb_basal_expr = symplified_model["symb_basal_expr"]
    init_cond = symplified_model["init_cond"]
    odes = symplified_model["odes"]
    lambdified_init_cond = lambdify_init_cond(param_to_estimate, init_cond)
    return Model(lambda t, p_VI, p_m1, p_m2, p_m4, p_m5, p_m6, p_HEb, p_Gb :
                    my_function(t,
                                param_to_estimate,
                                lambdified_init_cond,
                                odes,
                                **{"VI":p_VI, "m1":p_m1, "m2":p_m2, "m4":p_m4,
                                   "m5":p_m5, "m6":p_m6, "HEb":p_HEb, "Gb":p_Gb }))

# ## Parameter fitting for the liver compartment (ie. "Endogenous glucose production")
def generate_liver_model():
    '''
    Generates a model for the intestinal compartment.
    '''
    param_to_estimate = (kp1, kp2, kp3, kp4, ki, Gb)
    param_to_estimate_str = [str(p) for p in param_to_estimate]
    symplified_model = partial_evaluation(diabetic_values, param_to_estimate)
    params = symplified_model["parameter_values"]
    symb_basal_expr = symplified_model["symb_basal_expr"]
    init_cond = symplified_model["init_cond"]
    odes = symplified_model["odes"]
    lambdified_init_cond = lambdify_init_cond(param_to_estimate, init_cond)
    return Model(lambda t, p_kp1, p_kp2, p_kp3, p_kp4, p_ki, p_Gb :
                    my_function(t,
                                param_to_estimate,
                                lambdified_init_cond,
                                odes,
                                **{"kp1":p_kp1, "kp2":p_kp2, "kp3":p_kp3, "kp4":p_kp4,
                                   "ki":p_ki, "Gb":p_Gb }))

# ## Parameter fitting for the tissues compartment (ie. "utilization")
def generate_tissue_model():
    '''
    Generates a model for the intestinal compartment.
    '''
    param_to_estimate = (Fcns, Vm0, Vmx, Km0, p2U, Gb)
    param_to_estimate_str = [str(p) for p in param_to_estimate]
    symplified_model = partial_evaluation(diabetic_values, param_to_estimate)
    params = symplified_model["parameter_values"]
    symb_basal_expr = symplified_model["symb_basal_expr"]
    init_cond = symplified_model["init_cond"]
    odes = symplified_model["odes"]
    lambdified_init_cond = lambdify_init_cond(param_to_estimate, init_cond)
    return Model(lambda t, p_Fcns, p_Vm0, p_Vmx, p_Km0, p_p2U, p_Gb :
                    my_function(t,
                                param_to_estimate,
                                lambdified_init_cond,
                                odes,
                                **{"Fcns":p_Fcns, "Vm0":p_Vm0, "Vmx":p_Vmx, "Km0":p_Km0,
                                   "p2U": p_p2U, "Gb":p_Gb }))

# ## Parameter fitting for the pancreatic compartment (ie. "pancreatic secretion")
def generate_pancreatic_model():
    '''
    Generates a model for the intestinal compartment.
    '''
    param_to_estimate = (K, α, β, γ, Gb)
    param_to_estimate_str = [str(p) for p in param_to_estimate]
    symplified_model = partial_evaluation(diabetic_values, param_to_estimate)
    params = symplified_model["parameter_values"]
    symb_basal_expr = symplified_model["symb_basal_expr"]
    init_cond = symplified_model["init_cond"]
    odes = symplified_model["odes"]
    lambdified_init_cond = lambdify_init_cond(param_to_estimate, init_cond)
    return Model(lambda t, p_K, p_α, p_β, p_γ, p_Gb :
                    my_function(t,
                                param_to_estimate,
                                lambdified_init_cond,
                                odes,
                                **{"K":p_K, "α":p_α, "β":p_β, "γ":p_γ, "Gb":p_Gb}))

# ## Parameter fitting for the renal compartment (ie. "renal excretion")
def generate_renal_model():
    '''
    Generates a model for the intestinal compartment.
    '''
    param_to_estimate = (ke1, ke2, Gb)
    param_to_estimate_str = [str(p) for p in param_to_estimate]
    symplified_model = partial_evaluation(diabetic_values, param_to_estimate)
    params = symplified_model["parameter_values"]
    symb_basal_expr = symplified_model["symb_basal_expr"]
    init_cond = symplified_model["init_cond"]
    odes = symplified_model["odes"]
    lambdified_init_cond = lambdify_init_cond(param_to_estimate, init_cond)
    return Model(lambda t, p_ke1, p_ke2, p_Gb :
                    my_function(t,
                                param_to_estimate,
                                lambdified_init_cond,
                                odes,
                                **{"ke1":p_ke1, "ke2":p_ke2, "Gb":p_Gb}))

def Infer_seeker(p_tps, p_candidates, p_model, p_data, p_methods, p_compartment):
    '''
    Seeking for the best parameter for any model
    Inputs: - list of dates (p_tps)
            - list of parameters (p_candidates), lmfit object
            - model to fit (p_model), lmfit object
            - list of data (p_data)
            - list of methods for the inference (p_methods)
            - compartment name (p_compartment)
    Output: list for the best inference: - best lmfit output
                                         - best method
                                         - all lmfit output from all methods
    '''
    # result
    out = {}
    # computing the non-linear regression
    for mthd in p_methods:
        print(mthd)
        try:
            out_res = p_model.fit(p_data, p_candidates, t=p_tps, method=mthd)
            out[mthd] = out_res
            filestrparams = "results/%(compartment)s_%(name)s_params.p" % {"compartment": p_compartment, "name": mthd}
            filestrfit = "results/%(compartment)s_%(name)s_fit.p" % {"compartment": p_compartment, "name": mthd}
            filestrreport = "results/%(compartment)s_%(name)s_report.p" % {"compartment": p_compartment, "name": mthd}
            pickle.dump(out_res.best_values, open(filestrparams, 'wb'))
            pickle.dump(out_res.best_fit, open(filestrfit, 'wb'))
            pickle.dump(out_res.fit_report(), open(filestrreport, 'wb'))
        except ValueError:
            print(mthd+" is an invalid method for this session")
            pass
    return out

def plot_results(best_values, data):
    '''
    plots the results of a fitting against the data used to fit
    Inputs: - best_values: the attribute best_values of a modelResult
    - data: the simulated data used to fit.
    '''
    params = diabetic_values.copy()
    params.update({ eval(re.sub('p_', '', k)) : best_values[k]
                    for k in best_values.keys() })
    new_data = generate_simulated_data(params, nb_points=300)
    n = data[0].size
    figure_glucose = figure()
    plot(scipy.linspace(0, 300, 300), new_data[0], 'r-')
    plot(scipy.linspace(0, 300, n), data[0], 'r.')
    legend(["fitted glucose", "glucose data"])
    figure_insulin = figure()
    plot(scipy.linspace(0, 300, 300), new_data[1], 'b-')
    plot(scipy.linspace(0, 300, n), data[1], 'b.')
    legend(["fitted insulin", "plasma data"])

# methods = ['differential_evolution','basinhopping','ampgo','nelder','brute','leastsq','least_squares','lbfgsb','powell','cg','newton','cobyla','bfgs','tnc','trust-ncg','trust-exact','trust-krylov','trust-constr','dogleg','slsqp','emcee','shgo','dual_annealing']

# Workflow to fit data
# data = generate_simulated_data(healthy_values, nb_points=15)
# time_interv = scipy.linspace(0, 300, 15)
# %time results = Infer_seeker(time_interv, parameters, model_intestinal, data, ['leastsq'])

# to save the results:
# from lmfit.model import save_modelresult
# save_modelresult(results['leastsq'], 'leastsq_intestinal_results.sav')
# to load the results from a file:
# from lmfit.model import load_modelresult
# results = load_modelresult('leastsq_intestinal_results.sav')
