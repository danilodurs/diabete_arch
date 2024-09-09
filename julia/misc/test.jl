using BlackBoxOptim
include("../init/init.jl")
include("../param_estim.jl")
using .Dalla_Man_Init

fitted_vn = ["Gp","Ip","Ipo"]

time_points = Dict()
data = Dict()
sol = Dict()
for st in ["norm", "diab"]
    sol[st] = gen_solution(all_pn, all_p[st], all_vn, init_v[st], tspan, dalla_net!)
    (time_points[st], data[st]) = gen_data_set(30.0, all_vn, fitted_vn, tspan, sol[st])
end

(G,I,EGP,Ra,U,S) = (Dict(), Dict(), Dict(), Dict(), Dict(), Dict())
for st in ["norm", "diab"]
    (G[st],I[st],EGP[st],Ra[st],U[st],S[st]) = traces_of_obs_vars(sol[st],all_p[st])
end

compartments = Dict(
    "Glucose Kinetics"      => ["pk_1", "pk_2","pGb","pIb","pEGPb"],
    "Insulin Kinetics"      => ["pV_I", "pm_1", "pm_2", "pm_4", "pm_5", "pm_6", "pHE_b", "pGb","pIb","pEGPb"],
    "Ra"                    => ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc", "pα", "pβ","pGb","pIb","pEGPb"],
    "Ra2"                   => ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc", "pα", "pβ","pGb"],
    "Endogenous Production" => ["pk_p1", "pk_p2", "pk_p3", "pk_p4", "pk_i","pGb","pIb","pEGPb"],
    "Utilization"           => ["pF_cns", "pV_m0", "pV_mx", "pxK_m0", "pxK_mx", "pp_2U","pGb","pIb","pEGPb"],
    "Secretion"             => ["pK", "pxα_Y", "pxβ_Y", "pγ", "ph","pGb","pIb","pEGPb"],
    "Renal Excretion"       => ["pk_e1", "pk_e2","pGb","pIb","pEGPb"],
)

lower, upper = lower_upper_bounds(params, 1.5)

free_pn = compartments["Glucose Kinetics"]

prob = Dict()
res = Dict()

# création du ODEProblem
partial_model! = gen_model(free_pn, all_pn, all_p["diab"], dalla_net!)
params_GK = [all_p["diab"][free_pn[i]] for i in 1:length(free_pn)]
u0 = [ init_v["diab"][all_vn[i]] for i in 1:length(all_vn) ]

probGK = ODEProblem(partial_model!,u0,tspan,params_GK)

# estimation de paramètres
idx_of_free_pn = Dict(free_pn[i] => i for i in 1:length(free_pn))
free_p = [ all_p["diab"][pn] for pn in free_pn ]

function prob_update(prob,p)
    # This function is used to adjust the initial state of the ODE problem at
    # each iteration of the estimation algorithm. The initial state is
    # dependant on the (possibly estimated) parameters. 
    tmp = [(pn in free_pn) ? p[idx_of_free_pn[pn]] : all_p[pn]
            for pn in ["pF_cns", "pk_1", "pk_2", "pm_1", "pm_2", "pm_4", "pγ", "pHE_b", "pGb", "pIb", "pEGPb", "pV_G", "pV_I"] ]
    pF_cns, pk_1, pk_2, pm_1, pm_2, pm_4, pγ, pHE_b, pGb, pIb, pEGPb, pV_G, pV_I = tmp

    m30d = pHE_b*pm_1/(1-pHE_b)

    prob.u0[1]  = pGb*pV_G # G_p0
    prob.u0[2]  = (pF_cns-pEGPb+pk_1*prob.u0[1])/pk_2 # G_t0
    prob.u0[4]  = pIb*pV_I # I_p0
    prob.u0[3]  = prob.u0[4]*(pm_4+pm_2)/pm_1 # I_l0
    prob.u0[8]  = pIb # I_10
    prob.u0[9]  = pIb # I_d0
    prob.u0[12] = (prob.u0[4]*pm_4+prob.u0[3]*m30d)/pγ # I_po0

    return remake(prob,u0=convert.(eltype(p),prob.u0),p=p)
end


#####################################
# test de la somme des carrés sur avec l'inférence pour Ra et estimation basale
#####################################
using Serialization
include("./dalla_man.jl")

# using .Dalla_Man

include("./param_estim.jl")
include("./init.jl")

res = deserialize("results/Gp_with_basal_estimation/Ra results")

sub_pn = ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc", "pα", "pβ","pGb","pIb","pEGPb"]

# génération de la solution objective avec les paramètres "normaux"
init_v = gen_init_states(all_p["norm"])
sol_norm = gen_solution(all_pn, all_p["norm"], all_vn, init_v, tspan, dalla_net!)
(time_points, data) = gen_data_set(delay, all_vn, fitted_vn, tspan, sol_norm)

# génération de la solution objective avec les paramètres inférés
new_params = update_params(sub_pn, res.minimizer, all_pn, all_p["diab"])
init_v_inf = gen_init_states(new_params)
sol_inf = gen_solution(all_pn, new_params, all_vn, init_v_inf, tspan, dalla_net!)

# calcul du résiduel (somme des différences quadratiques) à partir de t=0
sumT0 = sum([(d-sol_inf(t)[1])^2 for (t,d) in zip(time_points,data[1,:])])
# cette somme doit aussi pouvoir s'écrire ainsi:
# sumT0 = sum(abs2, data[1,:]-sol_inf.(time_points)[1])

sumT1 = sumT0 - (data[1,1]-sol_inf(time_points[1])[1])^2

res.minimum # qui est pratiquement égal à sumT1 !!

#
# On recommence avec mais cette fois sans estimation basale
#

res = deserialize("results/Gp_without_basal_estimation/Ra results")

sub_pn = ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc", "pα", "pβ"]

# génération de la solution objective avec les paramètres "normaux"
init_v = gen_init_states(all_p["norm"])
sol_norm = gen_solution(all_pn, all_p["norm"], all_vn, init_v, tspan, dalla_net!)
(time_points, data) = gen_data_set(delay, all_vn, fitted_vn, tspan, sol_norm)

new_params = update_params(sub_pn, res.minimizer, all_pn, all_p["diab"])
init_v_inf = gen_init_states(new_params)
sol_inf = gen_solution(all_pn, new_params, all_vn, init_v_inf, tspan, dalla_net!)

sumT0 = sum([(d-sol_inf(t)[1])^2 for (t,d) in zip(time_points,data[1,:])])

sumT1 = sumT0 - (data[1,1]-sol_inf(time_points[1])[1])^2

res.minimum # qui est pratiquement égal à sumT1 !!

##############
# test de la somme des carrés sur l'exemple du tutoriel de DiffEqParamEstim
##############
using DifferentialEquations
using DiffEqParamEstim

function f(du,u,p,t)
    du[1] = dx = p[1]*u[1] - u[1]*u[2]
    du[2] = dy = -3*u[2] + u[1]*u[2]
  end
  
u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5]
prob = ODEProblem(f,u0,tspan,p)

sol = solve(prob,Tsit5())
t = collect(range(0,stop=10,length=200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)

cost_function = build_loss_objective(prob,Tsit5(),L2Loss(t,data),maxiters=10000,verbose=false)

lower = [0.0]
upper = [3.0]
result = optimize(cost_function, lower, upper, [0.42], Fminbox(BFGS()))

prob_inf = ODEProblem(f,u0,tspan,result.minimizer)
sol_inf = solve(prob_inf,Tsit5())

sum([(d-sol_inf(tp)[1])^2 for (tp,d) in zip(t,data[1,:])])+sum([(d-sol_inf(tp)[2])^2 for (tp,d) in zip(t,data[2,:])]) # quasiment égal à result.minimum

# idem mais en fittant seulement sur la première variable:
cost_function2 = build_loss_objective(prob,Tsit5(),L2Loss(t,data[1,:]),maxiters=10000,verbose=false,save_idxs=[1])
result2 = optimize(cost_function2, lower, upper, [1.42], Fminbox(BFGS()))