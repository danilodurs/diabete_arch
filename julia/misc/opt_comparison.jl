# Comparison of optimisation algorithms and packages

include("../models/dalla_man.jl")
include("../param_estim.jl")
include("../init/init.jl")

using .Init
using Serialization
using Plots
using NLopt
using Distributions
using BenchmarkTools
using BlackBoxOptim

model = dalla_net!

# compute init states for both normal and diabetic subject types
init_v = Dict(st => gen_init_states(all_p[st]) for st in ["norm","diab"])

# generate data points to fit 
sol = Dict(st => gen_solution(all_pn, all_p[st], all_vn, init_v[st], tspan, model) for st in ["norm","diab"])

time_points = Dict()
data = Dict()
for st in ["norm", "diab"]
    (time_points[st], data[st]) = gen_data_set(delay, all_vn, fitted_vn, tspan, sol[st])
    #(time_points[st], data[st]) = gen_data_set(delay, all_vn, all_vn, tspan, sol[st])
end

(G,I,EGP,Ra,U,S) = (Dict(), Dict(), Dict(), Dict(), Dict(), Dict())
for st in ["norm", "diab"]
    (G[st],I[st],EGP[st],Ra[st],U[st],S[st]) = traces_of_obs_vars(sol[st],all_p[st])
end

##########################
# parameter estimation
##########################
free_pn = ["pk_1", "pk_2", "pV_mx", "pxK_m0", "pp_2U"]

# first, fix common variables to all optimizations: lower, upper bounds, cost function, etc.
lower, upper = lower_upper_bounds(DM_params, coeff_bounds)
lower["pf"] = 0.2
upper["pf"] = 0.95
free_lower = [ lower[pn] for pn in free_pn ]
free_upper = [ upper[pn] for pn in free_pn ]

prob = Dict()
res = Dict()
figs = Dict()

idx_of_free_pn = Dict(free_pn[i] => i for i in 1:length(free_pn))

function prob_update(prob,p)
    # This function is used to adjust the initial state of the ODE problem at
    # each iteration of the estimation algorithm. The initial state is
    # dependant on the (possibly estimated) parameters. 
    tmp = [(pn in free_pn) ? p[idx_of_free_pn[pn]] : all_p["diab"][pn]
            for pn in ["pF_cns", "pk_1", "pk_2", "pm_1", "pm_5", "pm_6", "pγ", "pHE_b", "pGb", "pIb", "pEGPb", "pV_G", "pV_I"] ]
    pF_cns, pk_1, pk_2, pm_1, pm_5, pm_6, pγ, pHE_b, pGb, pIb, pEGPb, pV_G, pV_I = tmp
    m30d = pHE_b*pm_1/(1-pHE_b)
    pSb = Sb(pm_5,pm_6,pHE_b)
    pm_4 = m_4(pm_5,pm_6,pIb,pV_I,pHE_b)
    pm_2 = m_2(pm_5,pm_6,pIb,pV_I,pHE_b)
    # pSb = (pm_6-pHE_b)/pm_5
    # pm_4 = (2/5)*(pSb/(pIb*pV_I))*(1-pHE_b) # from eq (9)
    # pm_2 = (pSb/(pIb*pV_I)-(pm_4/(1-pHE_b)))*(1-pHE_b)/pHE_b # from eq (9)


    prob.u0[1]  = pGb*pV_G # G_p0
    prob.u0[2]  = (pF_cns-pEGPb+pk_1*prob.u0[1])/pk_2 # G_t0
    prob.u0[4]  = pIb*pV_I # I_p0
    prob.u0[3]  = prob.u0[4]*(pm_4+pm_2)/pm_1 # I_l0
    prob.u0[8]  = pIb # I_10
    prob.u0[9]  = pIb # I_d0
    #prob.u0[12] = (prob.u0[4]*pm_4+prob.u0[3]*m30d)/pγ # I_po0
    prob.u0[12] = pSb/pγ # I_po0

    return remake(prob,u0=convert.(eltype(p),prob.u0),p=p)
end


# generate cost function with initial ODEproblem instantiated with diabetic parameters and initial conditions
prob["initial"] = gen_prob(free_pn, all_pn, all_p["diab"], all_vn, init_v["diab"], tspan, model)
cost_function = build_loss_objective(prob["initial"],Tsit5(),L2Loss(time_points["norm"],data["norm"]),prob_generator = prob_update, maxiters=1000000,verbose = false, save_idxs = indexin(fitted_vn,all_vn))

free_p = [ all_p["diab"][pn] for pn in free_pn ] 

#####
# param estimation using Optim package with BFGS algorithm
free_p_estim = Dict()
losses = Dict()

@btime res["Optim_Fminbox(BFGS)"] = Optim.optimize(cost_function, free_lower, free_upper, free_p, Fminbox(BFGS()))
# 147s
free_p_estim["Optim_Fminbox(BFGS)"] = res["Optim_Fminbox(BFGS)"].minimizer
losses["Optim_Fminbox(BFGS)"] = res["Optim_Fminbox(BFGS)"].minimum

#####
# param estimation using Optim package with NelderMead algorithm
@btime res["Optim_Fminbox(NelderMead)"] = Optim.optimize(cost_function, free_lower, free_upper, free_p, Fminbox(NelderMead()))
# 45s
free_p_estim["Optim_Fminbox(NelderMead)"] = res["Optim_Fminbox(NelderMead)"].minimizer
losses["Optim_Fminbox(NelderMead)"] = res["Optim_Fminbox(NelderMead)"].minimum

#####
# param estimation using NLopt package with NelderMead algorithm
opt = Opt(:LN_NELDERMEAD, length(free_pn))
opt.min_objective = (p,grad) -> cost_function(p)
opt.xtol_rel = 0
opt.xtol_abs = 0
opt.lower_bounds = free_lower
opt.upper_bounds = free_upper

@btime res["NLopt_NelderMead"] = NLopt.optimize(opt, free_p)
# 579ms, has stopped on ":XTOL_REACHED" ??!!
free_p_estim["NLopt_NelderMead"] = res["NLopt_NelderMead"][2]
losses["NLopt_NelderMead"] = res["NLopt_NelderMead"][1]

# opt = Opt(:GN_DIRECT, length(free_pn))
# opt.min_objective = (p,grad) -> cost_function(p)
# opt.lower_bounds = free_lower
# opt.upper_bounds = free_upper

# @btime res["GN_DIRECT"] = NLopt.optimize(opt, free_p)
# trop long...

#####
# param estimation using BlackBoxOptim package
bounds = Tuple{Float64, Float64}[(l,u) for (l,u) in zip(free_lower,free_upper)]
@btime res["BlockBoxOptim"] = bboptimize(cost_function;SearchRange = bounds, MaxSteps = 11e3)
free_p_estim["BlockBoxOptim"] = best_candidate(res["BlockBoxOptim"])
losses["BlockBoxOptim"] = best_fitness(res["BlockBoxOptim"])
#####
# comparing curves

# plot normal and diabetic
plotly()

plot(sol["norm"].t,G["norm"],label="norm")
plot!(sol["diab"].t,G["diab"],label="diab")

for id in ["Optim_Fminbox(BFGS)", "Optim_Fminbox(NelderMead)", "NLopt_NelderMead", "BlockBoxOptim"]
    # all_p[id] = update_params(free_pn, free_p_estim[id], all_pn, all_p["diab"])
    # init_v[id] = gen_init_states(all_p[id])
    # sol[id] = gen_solution(all_pn, all_p[id], all_vn, init_v[id], tspan, model)
    # (G[id],I[id],EGP[id],Ra[id],U[id],S[id]) = traces_of_obs_vars(sol[id],all_p[id])
    plot!(sol[id].t,G[id],label=id)
end