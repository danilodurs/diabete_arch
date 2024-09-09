include("../models/dalla_man.jl")
include("../param_estim.jl")

using Serialization
using Plots
using Distributions
using LikelihoodProfiler

results_folder = "../results/PL"
model = dalla_net!
fitted_vn = ["Gp","Ip","Qgut"]
tspan = (0.0,420.0)
fitted_vn_idx = indexin(fitted_vn,all_vn)

# r = deserialize("./results/Gp_without_basal_estimation/Ra_results")

compartments = Dict(
     "GK"      => ["pk_1", "pk_2"],
     "IK"      => ["pV_I", "pm_1", "pm_5", "pm_6"], # we assume pHEb constant
     "Ra"      => ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc"],
     "EGP"     => ["pk_p2", "pk_p3", "pk_p4", "pk_i"],
     "U"       => ["pV_m0", "pV_mx", "pxK_m0", "pp_2U"], # we assume pF_cns and pxK_mx constant
     "S"       => ["pK", "pxα_Y", "pxβ_Y"], # we assume pγ constant
     "RE"      => ["pk_e1", "pk_e2"],
     "IS"      => ["pIb","pEGPb"] # init states
)

#####################################
# generate pseudo random data 
#####################################

# random data around normal data with 0.5% of standard deviation
pdev = 0.05
t = collect(range(tspan[1],stop=tspan[2],length=100))

# compute init states for both normal and diabetic subject types
init_v = Dict(st => gen_init_states(all_p[st]) for st in ["norm","diab"])

sol = Dict(st => gen_solution(all_pn, all_p[st], all_vn, init_v[st], tspan, model; saveat=t) for st in ["norm","diab"])

function generate_data(sol,t,fitted_vn_idx,pdev)
  randomized = VectorOfArray([[rand(Normal(sol(t[j])[i],max(abs(sol(t[j])[i])*pdev,0.5)),1)[1] for i in fitted_vn_idx] for j in 1:length(t)])
  convert(Array,randomized)
end

aggregate_data = Dict( st => convert(Array,VectorOfArray([generate_data(sol[st],t,fitted_vn_idx,pdev) for i in 1:200])) for st in ["norm","diab"])
distributions = Dict( st => [fit_mle(Normal,aggregate_data[st][i,j,:]) for i in 1:length(fitted_vn), j in 1:100] for st in ["norm","diab"])

#############################
# define cost function
#############################

# free_pn = [pn for pn in all_pn if all_p["norm"][pn] != all_p["diab"][pn] && pn != "pV_G" && pn != "pV_I"]
free_pn = Vector{String}()
for c in keys(compartments)
  append!(free_pn, compartments[c])
end

prob = Dict( st => gen_prob(free_pn,all_pn,all_p[st],all_vn,init_v[st],tspan,model) for st in ["norm", "diab"])

idx_of_free_pn = Dict(free_pn[i] => i for i in 1:length(free_pn))
free_p = Dict(st => [ all_p[st][pn] for pn in free_pn ] for st in ["norm", "diab"])

prob_updates = Dict( st => 
  function prob_update(prob,p)
    # This function is used to adjust the initial state of the ODE problem at
    # each iteration of the estimation algorithm. The initial state is
    # dependant on the (possibly estimated) parameters. 
    tmp = [(pn in free_pn) ? p[idx_of_free_pn[pn]] : all_p[st][pn]
      for pn in ["pF_cns", "pk_1", "pk_2", "pm_1", "pm_5", "pm_6", "pγ", "pHE_b", "pIb", "pEGPb", "pV_G", "pV_I", "pV_m0", "pxK_m0"] ]
    pF_cns, pk_1, pk_2, pm_1, pm_5, pm_6, pγ, pHE_b, pIb, pEGPb, pV_G, pV_I, pV_m0, pxK_m0 = tmp
    pGb = Gb(pk_1, pk_2, pV_m0, pxK_m0, pF_cns, pEGPb, pV_G)
    m30d = pHE_b*pm_1/(1-pHE_b)
    pSb = Sb(pm_5,pm_6,pHE_b)
    pm_4 = m_4(pm_5,pm_6,pIb,pV_I,pHE_b)
    pm_2 = m_2(pm_5,pm_6,pIb,pV_I,pHE_b)

    prob.u0[1]  = pGb*pV_G # G_p0
    prob.u0[2]  = (pF_cns-pEGPb+pk_1*prob.u0[1])/pk_2 # G_t0
    prob.u0[4]  = pIb*pV_I # I_p0
    prob.u0[3]  = prob.u0[4]*(pm_4+pm_2)/pm_1 # I_l0
    prob.u0[8]  = pIb # I_10
    prob.u0[9]  = pIb # I_d0
    prob.u0[12] = pSb/pγ # I_po0

    return remake(prob,u0=convert.(eltype(p),prob.u0),p=p)
  end
  for st in ["norm", "diab"] )

###############################
# define bounds and loss for profile likelihood
###############################
coeff_bounds = 2.0
theta_bounds = Dict()
scan_bounds = Dict()

for st in ["norm", "diab"]
  lower, upper = lower_upper_bounds(DM_params, coeff_bounds, st=st)
  lower["pf"] = 0.2
  upper["pf"] = 1.0

  # we need slightly wider intervals for theta bounds because scan bounds have to be narrower than theta bounds
  Δ_bounds = Dict(pn => abs(lower[pn]-upper[pn]) for pn in all_pn)
  wider_lower = Dict(pn => lower[pn]-Δ_bounds[pn]*1e-3 for pn in all_pn)
  wider_upper = Dict(pn => upper[pn]+Δ_bounds[pn]*1e-3 for pn in all_pn)

  # theta_bounds
  theta_bounds[st] = [ (wider_lower[pn],wider_upper[pn]) for pn in free_pn ]
  scan_bounds[st] = [ (lower[pn],upper[pn]) for pn in free_pn ]
end

interv = Dict(st => Vector{ParamInterval}(undef,length(free_p[st])) for st in ["norm", "diab"])

###################################################################
# computation of profile likelihoods w.r.t normal and diabetic data
###################################################################

for st in ["norm","diab"]
  # generate cost function on the basis of LogLikeLoss
  cost_function = build_loss_objective(prob[st],Tsit5(),LogLikeLoss(t,distributions[st]),prob_generator = prob_updates[st], maxiters=1000000,verbose = false, save_idxs = fitted_vn_idx)
  # log loss
  loss = cost_function(free_p[st])

  # log threshold
  α = loss + cquantile(Chisq(1), 0.05)
  # ou bien α = loss + cquantile(Chisq(1), 0.05) / 2 ??

  # compute confidence intervals
  for i in eachindex(free_p[st])
     interv[st][i] = get_interval(
        free_p[st],
        i,
        cost_function,
        :CICO_ONE_PASS;
        #:LIN_EXTRAPOL;
        loss_crit = α,
        theta_bounds = theta_bounds[st],
        scan_bounds = scan_bounds[st][i],
      )
  end
end

using DataFrames

resultats = Dict( st => DataFrame(
                          Parameters = free_pn,
                          StatusLower = [p.result[1].status for p in interv[st]],
                          StatusUpper = [p.result[2].status for p in interv[st]],
                          CILower = [p.result[1].value for p in interv[st]],
                          CIUpper = [p.result[2].value for p in interv[st]],
                          FittedValues = free_p[st]
                        ) 
                  for st in ["norm", "diab"] )

using CSV
for st in ["norm", "diab"]
  replace(resultats[st].StatusLower, nothing => missing)
  replace(resultats[st].StatusUpper, nothing => missing)
  resultats[st].CILower = replace(resultats[st].CILower, nothing => missing)
  resultats[st].CIUpper = replace(resultats[st].CIUpper, nothing => missing)
end
CSV.write("$results_folder/PL_results_norm.csv",resultats["norm"])
CSV.write("$results_folder/PL_results_diab.csv",resultats["diab"])

# iterpolation of profile likelihoods
for i in eachindex(free_p["norm"])
  update_profile_points!(interv["norm"][i])
  update_profile_points!(interv["diab"][i])
end

figs = Dict("norm" => Dict(), "diab" => Dict())
plotly()
for (pn,i) in zip(free_pn,1:length(free_pn))
  figs["norm"][pn] = plot(interv["norm"][i])
  savefig(figs["norm"][pn],"$results_folder/norm/$pn.png")
  figs["diab"][pn] = plot(interv["diab"][i])
  savefig(figs["diab"][pn],"$results_folder/diab/$pn.png")
end
