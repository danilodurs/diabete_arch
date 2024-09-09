using Catalyst, ModelingToolkit, IfElse, Latexify, DifferentialEquations, Plots, Serialization, CSV, StructuralIdentifiability
using CMAEvolutionStrategy, DiffEqParamEstim, LikelihoodProfiler, Distributions, DataFrames, StatsBase, ExcelFiles
include("../models/parameters.jl")
include("../myOptim/myCMAEvolutionStrategy.jl")
#include("../myOptim/myBlackBoxOptim.jl")
include("../misc/param_estim.jl")
cd("main")

## 

###############################
### data extraction           #
###############################


mpnorm  = load("../data/Minipigs_FullXylose.xlsx","Disparition D-xylose IV") |> DataFrame

timeticks = mpnorm[6:17,1]
idstr = Vector(mpnorm[3,2:end])
weights = Dict()
for i in 1:(length(idstr))
    weights[idstr[i]] = Meta.parse(mpnorm[4,(2+i-1)])
end
all_data_dx = Dict()
for i in 1:length(idstr)
    all_data_dx[idstr[i]] = mpnorm[6:end,(2+i-1)]
end 
for id in idstr
    Vg = all_prm[:DXDose][1]/(max(all_data_dx[id]...)*weights[id])
    all_data_dx[id] = Vg*all_data_dx[id]
end


# we add mean data row
tmp_matrix = Array{Float64,2}(undef,length(idstr),length(timeticks))
i = 1
for mp_id in idstr
    tmp_matrix[i,:] = all_data_dx[mp_id]
    i += 1
end
all_data_dx["MEAN"]=vec(mean(tmp_matrix; dims=1))
push!(idstr, "MEAN")

##

# import models

include("../models/elim_models.jl")

##

#######################
### MAIN
#######################

# 2 compartment model with elimination in remote compartment
u0 = [0.0,0.0]
results = run_all_estimations(crn_2comp_ke4,u0,idstr,timeticks,all_data_dx)
plts = make_plots(crn_2comp_ke4,u0,results,all_data_dx,timeticks)
save_plots(plts,"/home/lhoussai/tmp/figs/2comp_ke4_cmaes/")

##

# 2 compartment model without elimination in remote compartment
u0 = [0.0,0.0]
results = run_all_estimations(crn_2comp,u0,idstr,timeticks,all_data_dx)
plts = make_plots(crn_2comp,u0,results,all_data_dx,timeticks)
save_plots(plts,"/home/lhoussai/tmp/figs/2comp_ke4_cmaes/")

##

# 2 compartments model with 1st order elimination in remote compartment and no elimination in vascular compartment => THE BEST ONE
u0 = [0.0,0.0]
results = run_all_estimations(crn_2comp_remote_elim,u0,idstr,timeticks,all_data_dx)
plts = make_plots(crn_2comp_remote_elim,u0,results,all_data_dx,timeticks)
save_plots(plts,"/home/lhoussai/tmp/figs/2comp_remote_elim_cmaes/")

##
u0 = [0.0,0.0]
# let's remove the initial plateau from the data
all_data_dx_4 = Dict()
for minipig_id in idstr
    all_data_dx_4[minipig_id] = all_data_dx[minipig_id][4:end]
end 
timeticks_4 = timeticks[4:end] .- 15.0

results_sans_inj = run_all_estimations(crn_2comp_sans_inj, u0, idstr,timeticks_4,all_data_dx_4)
plts = make_plots(crn_2comp_sans_inj,u0,results_sans_inj,all_data_dx_4,timeticks_4)
save_plots(plts,"/home/lhoussai/tmp/figs/2comp_remote_sans_inj_cmaes/")

## 

# 1 compartment model with power exponetial elimination

# remove first 0 data point à time 0
all_data_dx_0 = Dict()
for minipig_id in idstr
    all_data_dx_0[minipig_id] = all_data_dx[minipig_id][2:end]
end 
timeticks_0 = timeticks[2:end] .- 5.0

results_ela = run_all_estimations(idstr,timeticks_0,all_data_dx_0)
plts = make_plots(results_ela,all_data_dx_0,timeticks_0)
save_plots(plts,"/home/lhoussai/tmp/figs/powexp_cmaes/")

## 

##
# Profile Likelihood

(cost_function, p0, lowerb, upperb) = prepare_estimation(crn_2comp_remote_elim,u0,timeticks,)

## 
# Structural Identifiability

ode_2comp_ke4 = convert(ODESystem, gen_crn_powexp(100))
measured_quantities = [y ~ DX_p]
local_id_all = assess_local_identifiability(ode_2comp_ke4, measured_quantities=measured_quantities, p=0.99)

glb_id_all = assess_identifiability(ode_2comp_ke4, measured_quantities=measured_quantities, p=0.99)

## 

##### Test
tspan= (0.0, timeticks_0[end])
crn = gen_crn_powexp(all_data_dx_0["AR7"][1])
prob = ODEProblem(convert(ODESystem,crn),[all_data_dx_0["AR7"][1]],tspan,zeros(length(parameters(crn))))


prob = remake(prob; p = results["AR7"][1], u0 = all_data_dx_0["AR7"][1])
sol = solve(prob)
    plt = plot(sol, xlabel="t", ylabel="mg/kg", label="DXp", vars=(0,1), title=fig_name*" (loss = "*string(round(loss))*")")
    plt = plot!(timeticks,data, m=:circle, label = "DX data")
    plt = annotate!(580, trunc(Int,maximum(data))-100, text("prm = "*string(round.(prm,sigdigits=3)), :grey, :right, 10))

plts[minipig_id] = plot_simulation(prob,results[minipig_id][1],results[minipig_id][2],data[minipig_id],timeticks,minipig_id)

##### TEST fin


##
## plot the result...

prob = remake(prob;p=mini)
new_sol = solve(new_prob)

plot(new_sol)
plot!(timeticks, all_data_dx[mn_id], m=:circle)

## Test powerexp model
ode_powerexp = convert(ODESystem,crn_powexp)
prob_powerexp = ODEProblem(ode_powerexp,[0.0],(0.1,200),[3.0,0.017])
sol_powerexp = solve(prob_powerexp)
plot(sol_powerexp)

k=0.02
β=0.5
ela(t) = 10*exp(-(k*t)^β)
plot!(ela, xlims=(0,100))