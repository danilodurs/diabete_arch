using Serialization
using Plots
using Distributions
using CSV
using DataFrames

include("../models/dalla_man.jl")
include("../param_estim.jl")
include("../init/init.jl")

using .Init


# method = "optim"
# include("../myOptim/myOptim.jl")
# include("../myOptim/myBlackBoxOptim.jl")
include("../myOptim/myCMAEvolutionStrategy.jl")

# results_folder = "../results/dataless/optim"
# results_folder = "../results/dataless/bboptim"
results_folder = "../results/dataless/cmaes"
model = dalla_net!

# compute init states for both normal and diabetic subject types
init_v = Dict(st => gen_init_states(all_p[st]) for st in ["norm","diab"])

# generate data points to fit 
sol = Dict(st => gen_solution(all_pn, all_p[st], all_vn, init_v[st], tspan, model) for st in ["norm","diab"])

time_points = Dict()
data = Dict()
for st in ["norm", "diab"]
    (time_points[st], data[st]) = gen_data_set(delay, all_vn, fitted_vn, tspan, sol[st])
end

(G,I,EGP,Ra,U,S) = (Dict(), Dict(), Dict(), Dict(), Dict(), Dict())
for st in ["norm", "diab"]
    (G[st],I[st],EGP[st],Ra[st],U[st],S[st]) = traces_of_obs_vars(sol[st],all_p[st])
end

# parameter estimation
lower, upper = lower_upper_bounds(DM_params, coeff_bounds)
lower["pf"] = 0.2
upper["pf"] = 0.95

prob = Dict()
res = Dict()

figs = Dict(cp => Dict() for cp in keys(compartments))

for comp in keys(compartments)
    prob[comp] = gen_prob(compartments[comp], all_pn, all_p["diab"], all_vn, init_v["diab"], tspan, model)
    res[comp] = param_estim(compartments[comp], all_p["diab"], fitted_vn, all_vn, lower, upper, time_points["norm"], data["norm"], prob[comp], my_optimize)
    
    all_p[comp] = update_params(compartments[comp], best_params(res[comp]), all_pn, all_p["diab"])
    init_v[comp] = gen_init_states(all_p[comp])
    sol[comp] = gen_solution(all_pn, all_p[comp], all_vn, init_v[comp], tspan, model)

    (G[comp],I[comp],EGP[comp],Ra[comp],U[comp],S[comp]) = traces_of_obs_vars(sol[comp],all_p[comp])

    for (var,svar,ylabel) in zip((G,I,EGP,Ra,U,S),("G","I","EGP","Ra","U","S"),("(mg/dL)", "(pmol/l)", "(mg/kg/min)", "(mg/kg/min)", "(mg/kg/min)", "(pmol/kg/min)"))
        figs[comp][svar] = plot(sol[comp].t, var[comp], title = svar, xlabel = "time (min)", ylabel = ylabel, labels = "Inferred", line = :dash, ylims = (0,maximum([maximum(var[comp]),maximum(var["diab"]),maximum(var["norm"])])))
        plot!(sol["diab"].t, var["diab"], labels = "T2D" )
        plot!(sol["norm"].t, var["norm"], labels = "Normal")
    end

    figs[comp]["all"] = plot(figs[comp]["G"],figs[comp]["I"],figs[comp]["EGP"],figs[comp]["Ra"],figs[comp]["U"],figs[comp]["S"],size=(1000,800) ,layout = (3,2))
    savefig(figs[comp]["all"],results_folder*"/$comp.png")
end

# save infered parameter values for each compartment
for cn in keys(compartments)
    lv = loss_value(res[cn]) 
    res_tmp = DataFrame(Param = compartments[cn], vals = best_params(res[cn])) 
    push!(res_tmp, ["loss", lv])
    CSV.write(results_folder*"/$cn results.csv",res_tmp)
    # serialize("results/$cn results",res[cn])
end

# generate bar plot
compartments_names = [cn for cn in keys(compartments)]
results = [loss_value(res[cn]) for cn in compartments_names]
comparaisons = bar(compartments_names,results,labels = nothing ,title = "L2 Minimization by Biological Function" )
savefig(comparaisons,results_folder*"/l2loss.png")
