from model import *
from parametersValues import parameter_values
from customSubstitutions import *

# The steady state object is made of
# - "conditions" used to substitute some symbols to actual values stating the steady state
# - "equations" which are the equations of the original model to which the conditions are applied
# - "symbolic_basal_values" the symbolic basal values (i.e. expressed in terms of symbols)
# - "numerical_basal_values" the numeric basal values (i.e. the evaluation of symbolic basal
#    values for some given parameter set et initial conditions)

def init_steady_state_analysis():
    """
        Initializes the steady state object.
        
        Returns the initialized steady state object.
    """
    steady_state = EmptyObject()
    steady_state.conditions = EmptyObject()
    steady_state.equations = model.equations.copy()
    steady_state.symbolic_basal_values = dict()
    steady_state.variables_to_basal = { x:y for (x,y) in list(zip(symbols.variables, symbols.variables_basal)) }
    return steady_state

def generate_steady_state_eqns(steady_state):
    """
        Generates the conditions for steady state and the subsequent equations.

        Let us first state the conditions of steady state:
        - no stimuli (D=0) and
        - each variable's derivative is set to zero
        - each variable is set to its basal value

        The ss equations are obtained by applying the conditions of steady
        state to the original model's equations.
        In addition, we state that h is equal to Gb as in Dalla Mann's paper.
    """
    steady_state.conditions = { D: 0 }
    for ev in model.symbols.variables:
        steady_state.conditions.update({ ev.diff(): 0 })
    steady_state.conditions.update(steady_state.variables_to_basal)
    steady_state.equations = subs(steady_state.equations, steady_state.conditions)
    steady_state.equations.update({ h : sympy.Eq(h, Gb) })

################################################################################
# Some usefull fonction to automate the resolution of the steady state equations
################################################################################

def solve_partial_equations(steady_state, variables):
    """
        Solves the equations eqns for the given variables.

        Returns a mapping from basal variables to expressions.
    """
    tmp_equations = [
        steady_state.equations[e]
        for e in variables
        if e in steady_state.equations.keys()
    ]
    return sympy.solve(tmp_equations, {steady_state.variables_to_basal[e] for e in variables})


def partial_update_of_steady_state(steady_state, sol):
    """
        Updates the symbolic basal values with new expressions given by the sol parameter:
        If symbolic_basal_values = { X->E for some X }
        returns the composition: sol o { X->sol(E) for the same X }
    """
    for k in list(steady_state.symbolic_basal_values.keys()):
        if type(steady_state.symbolic_basal_values[k])!=int:
            steady_state.symbolic_basal_values.update(
                {k:subs(steady_state.symbolic_basal_values[k],sol)}
            )
    steady_state.symbolic_basal_values.update(sol)

def partial_update_of_ss_equations(steady_state):
    """
        Returns the equations eqns updated with symbolic_basal_values.
    """
    return subs(steady_state.equations,steady_state.symbolic_basal_values)

def generate_symbolic_basal_values(steady_state):
    """
        Solves the steady state equations by symbolic manipulation.
    """
    # We can safely assume that that glucose renal excretion E is null at basal steady state.
    # This, by the way, brakes down the piecewise definition of the equation for Eb.
    steady_state.symbolic_basal_values.update({Eb:0})
    steady_state.equations.update({E : sympy.Eq(Eb, 0)})

    steady_state.equations = subs(steady_state.equations, {Eb : 0, h: Gb})

    # We first solve the trivial equations for G, HE, I, kempt, Qgut, Qsto, Qsto1, Qsto2, Ra, S, Uii, X, Y, m3 and EGP

    tmp_solutions = solve_partial_equations(
                        steady_state,
                        [Gp,HE,I,kempt,Qgut,Qsto,Qsto1,Qsto2,Ra,S,Uii,X,Y,m3,EGP]
                    )
    partial_update_of_steady_state(steady_state, tmp_solutions[0])
    steady_state.equations = partial_update_of_ss_equations(steady_state)

    # Solve sub-system for I1b, Idb and Ipb

    tmp_solutions = solve_partial_equations(steady_state,[I1,Id,Ip])
    partial_update_of_steady_state(steady_state, tmp_solutions)
    steady_state.equations = partial_update_of_ss_equations(steady_state)

    # Solve sub-system for Il, Ipo and Spo

    tmp_solutions = solve_partial_equations(steady_state,[Il,Ipo,Spo])
    partial_update_of_steady_state(steady_state, tmp_solutions[0])
    steady_state.equations = partial_update_of_ss_equations(steady_state)

    # Solve the sub-system for Ub and Uidb

    tmp_solutions = solve_partial_equations(steady_state,[U,Uid])
    partial_update_of_steady_state(steady_state, tmp_solutions)
    steady_state.equations = partial_update_of_ss_equations(steady_state)
    
    # Solve the sub-system for Spob
    # CAUTION: two solutions are found for Spob. For the parameter values given
    # in DallaMan paper only the 1st one gives a positive value for Gtb.
    # In principle, we should keep all possible symbolic values. To do later...
    
    tmp_solutions = sympy.solve(steady_state.equations[G],{Spob})
    partial_update_of_steady_state(steady_state, {Spob:tmp_solutions[0]})
    steady_state.equations = partial_update_of_ss_equations(steady_state)
    
    # Solve the sub-system for Gtb
    # CAUTION: again, here we have two possible symbolic solutions, we keep 
    # the 2nd one. To do: generalize...
    
    tmp_solutions = sympy.solve(steady_state.equations[Gt],{Gtb})
    partial_update_of_steady_state(steady_state, {Gtb:tmp_solutions[1]})
    steady_state.equations = partial_update_of_ss_equations(steady_state)