using DataFrames
using Serialization
using ExcelFiles
using Statistics
using Distributions
using CSV
using Plots

include("../models/dalla_man.jl")
include("../param_estim.jl")

#############################################
# Various initializations
#############################################

# optimization method

include("../myOptim/myOptim.jl")
# include("../myOptim/myBlackBoxOptim.jl")
# include("../myOptim/myCMAEvolutionStrategy.jl")

results_folder = "../results/obediab/optim"
# results_folder = "../results/obediab/bboptim"
# results_folder = "../results/obediab/cmaes"

dallaman_diab = "T2D model"
dallaman_norm = "normal model"


compartments = Dict(
     "GK"      => ["pk_1", "pk_2"],
     "IK"      => ["pV_I", "pm_1", "pm_5", "pm_6"], # we assume pHEb constant
     "Ra"      => ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc"],
     "EGP"     => ["pk_p2", "pk_p3", "pk_p4", "pk_i"],
     "U"       => ["pV_m0", "pV_mx", "pxK_m0", "pp_2U"], # we assume pF_cns and pxK_mx constant
     "S"       => ["pK", "pxα_Y", "pxβ_Y"], # we assume pγ constant
     "RE"      => ["pk_e1", "pk_e2"],
     "IS"      => ["pEGPb", "pIb"] #, "pGb","pIb"] # init states
)

fitted_vn = ["Gp","Ip"]

free_pn = Vector{String}()
for c in keys(compartments)
  append!(free_pn, compartments[c])
end

tspan = (0.0,420.0)


# the glucose of OBEDIAB is about 22g of glucose and 32g of amidon
DM_params["pD"] = [54000.0,54000.0]

distdA = deserialize("../data/distdA")
distdB = deserialize("../data/distdB")

dataA_mean = deserialize("../data/dataA_mean")
dataB_mean = deserialize("../data/dataB_mean")
ob_d_A_Gp = dataA_mean[1,:]
ob_d_A_Ip = dataA_mean[2,:]
ob_d_B_Gp = dataB_mean[1,:]
ob_d_B_Ip = dataB_mean[2,:]

BWs_A = deserialize("../data/weightsA")
BW_A = mean([BWs_A[id][1] for id in keys(BWs_A)])

BWs_B = deserialize("../data/weightsB")
BW_B = mean([BWs_B[id][1] for id in keys(BWs_B)])

################################################################
# Simulate data with normal and diabetic paramaters of Dalla Man
################################################################

t= [0.0,30.0,60.0,90.0,120.0,180.0]
time_points = Dict()
time_points["diab"] = t
time_points["norm"] = t

# compute init states for both normal and diabetic subject types
init_v = Dict(st => gen_init_states(all_p[st]) for st in ["norm","diab"])

# generate data points to fit 
sol = Dict(st => gen_solution(all_pn, all_p[st], all_vn, init_v[st], tspan, dalla_net!) for st in ["norm","diab"])

(G,I,EGP,Ra,U,S) = (Dict(), Dict(), Dict(), Dict(), Dict(), Dict())
for st in ["norm", "diab"]
    (G[st],I[st],EGP[st],Ra[st],U[st],S[st]) = traces_of_obs_vars(sol[st],all_p[st])
end

# parameter estimation
DM_params["BW"] = [BW_A,BW_A]
lower, upper = lower_upper_bounds(DM_params, 1.5)
lower["pf"] = 0.5
upper["pf"] = 1.0
lower["pk_abs"] = 0.01 # otherwise 0 is computed

all_p["start"] = Dict(pn => all_p["diab"][pn] for pn in all_pn)

#all_p["start"]["pGb"] = ob_d_A_Gp[1]/all_p["start"]["pV_G"]
all_p["start"]["pIb"] = ob_d_A_Ip[1]/all_p["start"]["pV_I"]

init_v["start"] = gen_init_states(all_p["start"])

dataA = Array{Float64,2}(hcat(ob_d_A_Gp, ob_d_A_Ip)')
prob = gen_prob(free_pn, all_pn, all_p["start"], all_vn, init_v["start"], tspan, dalla_net!)
resA = param_estim(free_pn, all_p["start"], fitted_vn, all_vn, lower, upper, time_points["diab"], distdA, prob, my_optimize)

# plotting results
comp = "visitA"

figs = Dict()
figs[comp] = Dict()

estimated_params = best_params(resA)
all_p[comp] = update_params(free_pn, estimated_params, all_pn, all_p["start"])

init_v[comp] = gen_init_states(all_p[comp])
sol[comp] = gen_solution(all_pn, all_p[comp], all_vn, init_v[comp], tspan, dalla_net!)

(G[comp],I[comp],EGP[comp],Ra[comp],U[comp],S[comp]) = traces_of_obs_vars(sol[comp],all_p[comp])

VarGly3A = ob_d_A_Gp / DM_params["pV_G"][2]
VarIns3A = ob_d_A_Ip / DM_params["pV_I"][2]

for (var,svar,ylabel) in zip((G,I,EGP,Ra,U,S),("G","I","EGP","Ra","U","S"),("(mg/dL)", "(pmol/l)", "(mg/kg/min)", "(mg/kg/min)", "(mg/kg/min)", "(pmol/kg/min)"))
    figs[comp][svar] = plot(sol[comp].t, var[comp], title = svar, linecolor = :blue, xlabel = "time (min)", ylabel = ylabel, labels = "\"Visit A\" model", ylims = (0,maximum([maximum(var[comp]),maximum(var["diab"])])))
    if(svar == "G")
       plot!(t,VarGly3A, linecolor = :red, labels = "\"Visit A\" data")
    #    plot!(sol["diab"].t, var["diab"], labels = "T2D model",ylims = (0,maximum([maximum(var[comp]),maximum(var["diab"]),maximum(VarGly3A)])))
       plot!(sol["diab"].t, var["diab"], labels = dallaman_diab, linecolor = :black, ylims = (0,maximum([maximum(var[comp]),maximum(var["diab"]),maximum(VarGly3A)])))
    elseif(svar == "I")
        plot!(t,VarIns3A, linecolor = :red, labels = "\"Visit A\" data")
        # plot!(sol["diab"].t, var["diab"], labels = "T2D model",ylims = (0,maximum([maximum(var[comp]),maximum(var["diab"]),maximum(VarIns3A)])))
        plot!(sol["diab"].t, var["diab"], labels = dallaman_diab, linecolor = :black, ylims = (0,maximum([maximum(var[comp]),maximum(var["diab"]),maximum(VarIns3A)])))
    else
        # plot!(sol["diab"].t, var["diab"], labels = "T2D model")
        plot!(sol["diab"].t, var["diab"], labels = dallaman_diab, linecolor = :black)
    end
end

figs[comp]["all"] = plot(figs[comp]["G"],figs[comp]["I"],figs[comp]["EGP"],figs[comp]["Ra"],figs[comp]["U"],figs[comp]["S"],size=(1000,800) ,layout = (3,2))
savefig(figs[comp]["all"],results_folder*"/$(comp).png")

lv = loss_value(resA) 
resultsA = DataFrame(Param = free_pn, vals = estimated_params) 
push!(resultsA, ["loss", lv])

CSV.write(results_folder*"/resultsA.csv",resultsA)

###############################################
# Now we estimate all the parameters to match the visit B starting from the parameters we just estimated from visit A
###############################################
DM_params["BW"] = [BW_B,BW_B]

all_p["start"] = copy(all_p["visitA"])
all_p["start"]["pIb"] = ob_d_B_Ip[1]/all_p["start"]["pV_I"]

init_v["start"] = gen_init_states(all_p["start"])

dataB = Array{Float64,2}(hcat(ob_d_B_Gp, ob_d_B_Ip)')

# for comp in keys(compartments)
prob = gen_prob(free_pn, all_pn, all_p["start"], all_vn, init_v["start"], tspan, dalla_net!)
resB = param_estim(free_pn, all_p["start"], fitted_vn, all_vn, lower, upper, time_points["diab"], distdB, prob, my_optimize)

# plotting results
comp = "visitB"

figs[comp] = Dict()

estimated_params = best_params(resB)
all_p[comp] = update_params(free_pn, estimated_params, all_pn, all_p["start"])

init_v[comp] = gen_init_states(all_p[comp])
sol[comp] = gen_solution(all_pn, all_p[comp], all_vn, init_v[comp], tspan, dalla_net!)

(G[comp],I[comp],EGP[comp],Ra[comp],U[comp],S[comp]) = traces_of_obs_vars(sol[comp],all_p[comp])

VarGly3B = ob_d_B_Gp / DM_params["pV_G"][2]
VarIns3B = ob_d_B_Ip / DM_params["pV_I"][2]

for (var,svar,ylabel) in zip((G,I,EGP,Ra,U,S),("G","I","EGP","Ra","U","S"),("(mg/dL)", "(pmol/l)", "(mg/kg/min)", "(mg/kg/min)", "(mg/kg/min)", "(pmol/kg/min)"))
    figs[comp][svar] = plot(sol[comp].t, var[comp], title = svar, linecolor=:green, xlabel = "time (min)", ylabel = "G (mg/dL)", labels = "\"Visit B\" model", ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"])])))
    if(svar == "G")
        plot!(t,VarGly3B, linecolor=:red, labels = "\"Visit B\" data")
        plot!(sol["visitA"].t, var["visitA"], linecolor=:blue, labels = "\"Visit A\" model",ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"]),maximum(VarGly3B)])))
    elseif(svar == "I")
        plot!(t,VarIns3B, linecolor=:red, labels = "\"Visit B\" data")
        plot!(sol["visitA"].t, var["visitA"], linecolor=:blue, labels = "\"Visit A\" model",ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"]),maximum(VarIns3B)])))
    else
        plot!(sol["visitA"].t, var["visitA"], linecolor=:blue, labels = "\"Visit A\" model")
    end
end

figs[comp]["all"] = plot(figs[comp]["G"],figs[comp]["I"],figs[comp]["EGP"],figs[comp]["Ra"],figs[comp]["U"],figs[comp]["S"],size=(1000,800) ,layout = (3,2))
savefig(figs[comp]["all"],results_folder*"/$(comp).png")

lv = loss_value(resB) 
resultsB = DataFrame(Param = free_pn, vals = estimated_params) 
push!(resultsB, ["loss", lv])

CSV.write(results_folder*"/resultsB.csv",resultsB)

###############################################
# parameter estimation compartment by compartment from estimated visitA to match visitB
###############################################

compartments = Dict(
     "GK+U" => ["pk_1", "pk_2", "pV_m0", "pV_mx", "pxK_m0", "pp_2U", "pIb", "pEGPb"],
     "IK+S" => ["pV_I", "pm_1", "pm_5", "pm_6", "pK", "pxα_Y", "pxβ_Y", "pIb", "pEGPb"], # we assume pHEb constant
     "Ra"   => ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc", "pIb", "pEGPb"],
     "EGP"  => ["pk_p2", "pk_p3", "pk_p4", "pk_i", "pIb", "pEGPb"],
     "RE"   => ["pk_e1", "pk_e2", "pIb", "pEGPb"],
)

for comp in keys(compartments)
    figs[comp] = Dict()
end

probs = Dict()
res = Dict()

dataB = Array{Float64,2}(hcat(ob_d_B_Gp, ob_d_B_Ip)')

for comp in keys(compartments)
    probs[comp] = gen_prob(compartments[comp], all_pn, all_p["visitA"], all_vn, init_v["visitA"], tspan, dalla_net!)
    res[comp] = param_estim(compartments[comp], all_p["visitA"], fitted_vn, all_vn, lower, upper, time_points["diab"], distdB, probs[comp], my_optimize)

    estimated_params = best_params(res[comp])
    all_p[comp] = update_params(compartments[comp], estimated_params, all_pn, all_p["visitA"])
    init_v[comp] = gen_init_states(all_p[comp])
    sol[comp] = gen_solution(all_pn, all_p[comp], all_vn, init_v[comp], tspan, dalla_net!)

    (G[comp],I[comp],EGP[comp],Ra[comp],U[comp],S[comp]) = traces_of_obs_vars(sol[comp],all_p[comp])

    VarGly3 = ob_d_B_Gp / DM_params["pV_G"][2]
    VarIns3 = ob_d_B_Ip / DM_params["pV_I"][2]

    for (var,svar,ylabel) in zip((G,I,EGP,Ra,U,S),("G","I","EGP","Ra","U","S"),("(mg/dL)", "(pmol/l)", "(mg/kg/min)", "(mg/kg/min)", "(mg/kg/min)", "(pmol/kg/min)"))
        # figs[comp][svar] = plot(sol[comp].t, var[comp], linecolor = :green, title = svar, xlabel = "time (min)", ylabel = ylabel, labels = "\"Visit B\" model", ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"])])))
        figs[comp][svar] = plot(sol[comp].t, var[comp], linecolor = :blue, linestyle = :dash, title = svar, xlabel = "time (min)", ylabel = ylabel, labels = "\"Visit B\" model", ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"])])))
        if(svar == "G")
            # plot!(sol["visitA"].t, var["visitA"], linecolor = :blue, labels = "\"Visit A\" model",ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"]),maximum(VarGly3)])))
            plot!(sol["visitA"].t, var["visitA"], linecolor = :red, labels = "\"Visit A\" model",ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"]),maximum(VarGly3)])))
            # plot!(t,VarGly3, linecolor = :red, labels = "\"Visit B\" data")
            plot!(t,VarGly3, linecolor = :green, labels = "\"Visit B\" data")
         elseif(svar == "I")
            # plot!(sol["visitA"].t, var["visitA"], linecolor = :blue, labels = "\"Visit A\" model",ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"]),maximum(VarIns3)])))
            plot!(sol["visitA"].t, var["visitA"], linecolor = :red, labels = "\"Visit A\" model",ylims = (0,maximum([maximum(var[comp]),maximum(var["visitA"]),maximum(VarIns3)])))
            # plot!(t,VarIns3, linecolor = :red, labels = "\"Visit B\" data")
            plot!(t,VarIns3, linecolor = :green, labels = "\"Visit B\" data")
         else
            #  plot!(sol["visitA"].t, var["visitA"], linecolor = :blue, labels = "\"Visit A\" model")
             plot!(sol["visitA"].t, var["visitA"], linecolor = :red, labels = "\"Visit A\" model")
         end     
    end    

    figs[comp]["all"] = plot(figs[comp]["G"],figs[comp]["I"],figs[comp]["EGP"],figs[comp]["Ra"],figs[comp]["U"],figs[comp]["S"],size=(1000,800) ,layout = (3,2))
    savefig(figs[comp]["all"],results_folder*"/$comp.png")
end

# save parameter values in csv file
for comp in keys(compartments)
    local estimated_params = best_params(res[comp])
    local lv = loss_value(res[comp])
    DFresults = DataFrame(Param = compartments[comp], vals = estimated_params)
    push!(DFresults, ["loss", lv])

    CSV.write(results_folder*"/$comp.csv",DFresults)
end

# save loss results
compartments_names = [n for n in keys(compartments)]
results = [loss_value(res[comp]) for comp in compartments_names]
comparaisons = bar(compartments_names,results,labels = nothing ,title = "Log Likelihood Minimization by Biological Function" )
savefig(comparaisons,results_folder*"/loss.png")


f_iks = plot(figs["IK+S"]["G"],figs["IK+S"]["I"],figs["IK+S"]["Ra"],size=(1100,150) ,layout = (1,3))
f_ra = plot(figs["Ra"]["G"],figs["Ra"]["I"],figs["Ra"]["Ra"],size=(1100,150) ,layout = (1,3))
f_visita = plot(figs["visitA"]["G"],figs["visitA"]["I"],figs["visitA"]["Ra"],size=(1100,150) ,layout = (1,3))
f_visitb = plot(figs["visitB"]["G"],figs["visitB"]["I"],figs["visitB"]["Ra"],size=(1100,150) ,layout = (1,3))

savefig(f_iks,results_folder*"/IKS_GIRa.png")
savefig(f_ra,results_folder*"/Ra_GIRa.png")
savefig(f_visita,results_folder*"/VisitA_GIRa.png")
savefig(f_visitb,results_folder*"/VisitB_GIRa.png")

