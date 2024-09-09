using DifferentialEquations
using DiffEqParamEstim

"""
    gen_model(free_pn::Vector{String}, all_pn::Vector{String}, all_p::Dict{String,Float64}, model!::Function)::Function

given a model, this function generates a (partial) model for some free parameter names listed in `free_pn` 
(resp. `all_pn`), all parameter name/value maps of diabetic subjects collected 
in the dictionary `all_p` and basic functional model `model!`.
"""
function gen_model(free_pn::Vector{String},     # free parameter names
                   all_pn::Vector{String},      # all parameter names
                   all_p::Dict{String,Float64}, # parameter -> value dict
                   model!                       # model with all parameter free
                  )::Function
    idx_of_free_pn = indexin(free_pn,all_pn)
    # pp is the vector of all fixed param values expect the free parameters that are updated with p in partial_model!
    pp = [ all_p[pn] for pn in all_pn ]
    function partial_model!(du,u,p,t)
        for i in 1:length(idx_of_free_pn)
            pp[idx_of_free_pn[i]] = p[i]
        end
        model!(du,u,pp,t)
    end
    partial_model!
end # gen_model

"""
    gen_solution(all_pn::Vector{String}, all_p::Dict{String,Float64},
    all_vn::Vector{String},init_v::Dict{String,Float64},tspan::Tuple{Float64,Float64},
    model!::Function)::OrdinaryDiffEq.ODECompositeSolution

generates a solution to the ODE problem generated from given model, parameters and initial states.
"""
function gen_solution(
                      all_pn::Vector{String},       # parameter names: gives the order of parameters
                      all_p::Dict{String,Float64},  # parameter -> value dict
                      all_vn::Vector{String},       # variable names: gives the order of variables
                      init_v::Dict{String,Float64}, # variable -> init state dict
                      tspan::Tuple{Float64,Float64},# span ot time simulation
                      model!::Function;              # model with all parameter free
                      kwargs...
                     )::OrdinaryDiffEq.ODECompositeSolution
    params = [all_p[pn] for pn in all_pn]
    u0     = [init_v[vn] for vn in all_vn]
    prob   = ODEProblem(model!,u0,tspan,params)
    solve(prob;kwargs...)
end # gen_solution

"""
    gen_data_set(step::Float64,all_vn::Vector{String},tspan::Tuple{Float64,Float64},sol)::Tuple{Array{Float64,1},Array{Float64,2}}

generates the time course data set of the variables Gp, Ip and Ipo for a given model, 
parameters and initial states. Data points are given each `st` time step.
"""
function gen_data_set(
                      step::Float64,                # time step 
                      all_vn::Vector{String},       # variable names: gives the order of variables
                      fitted_vn::Vector{String},    # Variables which we want to extract data
                      tspan::Tuple{Float64,Float64},# span ot time simulation
                      sol::OrdinaryDiffEq.ODECompositeSolution # solution of an ODE problem
                      )::Tuple{Vector{Float64},Array{Float64,2}}
    
    idx_of_fitted_vn = indexin(fitted_vn, all_vn)
    time_points = collect(tspan[1]:step:tspan[2])
    data = [[sol(t)[idx] for t in time_points] for idx in idx_of_fitted_vn]
    data = convert(Array, transpose(hcat(data...)))
    return (time_points, data)
end # gen_data_set

"""
    gen_prob(free_pn::Array{String}, all_pn::Vector{String}, all_p::Dict{String,Float64}, 
    all_vn::Vector{String}, init_v::Dict{String,Float64}, tspan::Tuple{Float64,Float64}, 
    model!::Function)::ODEProblem

generates a problem for free parameters names listed in `free_pn`. By default, non free parameter values are fixed with diabetic values.
"""
function gen_prob(free_pn::Vector{String},      # free parameter names
                  all_pn::Vector{String},       # parameter names: gives the order of parameters
                  all_p::Dict{String,Float64},  # parameter -> value dict
                  all_vn::Vector{String},       # variable names: gives the order of variables
                  init_v::Dict{String,Float64}, # variable -> init state dict
                  tspan::Tuple{Float64,Float64}, # span ot time simulation
                  model!::Function              # model with all parameter free
                 )::ODEProblem
    
                 # problem's arguments
    partial_model!::Function = gen_model(free_pn,all_pn,all_p,model!)
    params::Array{Float64}  = [ all_p[free_pn[i]] for i in 1:length(free_pn) ]
    u0::Array{Float64}      = [ init_v[all_vn[i]] for i in 1:length(all_vn) ]

    return ODEProblem(partial_model!,u0,tspan,params)
end # gen_prob

"""
    lower_upper_bounds(params::Dict{String,Vector{Float64}}, # maps each parameter name to pairs of parameters values
                       coef::Float64;                             # determines the width of upper and lower bound
                       st = "diab"                           # specifies around which subject type parameters the interval have to be built
                      )::Tuple{Vector{Float64},Vector{Float64}}
"""
function lower_upper_bounds(params::Dict{String,Vector{Float64}}, # maps each parameter name to pairs of parameters values
                            coef::Float64;                        # determines the width of upper and lower bound
                            st = "diab"                           # specifies around which subject type parameters the interval have to be built
                           )::Tuple{Dict{String,Float64},Dict{String,Float64}}
    params_Δ = Dict(pn => coef*abs(params[pn][2] - params[pn][1]) for pn in keys(params))
    if (st == "diab") # then lower and upper bounds are given around diabetic parameters
        lower::Dict{String,Float64} = Dict(pn => max(0,params[pn][2]-params_Δ[pn]) for pn in keys(params))
        upper::Dict{String,Float64} = Dict(pn => params[pn][2] + params_Δ[pn] for pn in keys(params))
        # if normal parameters are not reachable we raise an exception
        for pn in keys(params)
            if !(lower[pn] <= params[pn][1] <= upper[pn])
                throw(DomainError(pn,"normal parameters are not reachable!"))
            end
        end
    else # otherwise they are given around normal parameters
        lower = Dict(pn => max(0,params[pn][1]-params_Δ[pn]) for pn in keys(params))
        upper = Dict(pn => params[pn][1] + params_Δ[pn] for pn in keys(params))
        # if diabetic parameters are not reachable we raise an exception
        for pn in keys(params)
            if !(lower[pn] <= params[pn][2] <= upper[pn])
                throw(DomainError(pn,"normal parameters are not reachable!"))
            end
        end
    end
    (lower,upper)
end # lower_upper_bounds

"""
    get_prob_update(free_pn, # free parameter names
                    all_p)   # parameter -> value dict
    
    return the problem update function used at each iteration of the parameter estimation.
"""
function get_prob_update(free_pn, all_p, all_vn)
    idx_of_free_pn = Dict(free_pn[i] => i for i in 1:length(free_pn))
    all_p = copy(all_p)

    function prob_update(prob,p)
        #the_pns = ["pF_cns", "pk_1", "pk_2", "pm_1", "pm_5", "pm_6", "pγ", "pHE_b", "pIb", "pEGPb", "pV_G", "pV_I", "pV_m0", "pxK_m0", "pD"]
        for i in 1:length(p)
            all_p[free_pn[i]] = p[i]
        end
        #params_dict = Dict(pn => (pn in free_pn) ? p[idx_of_free_pn[pn]] : all_p[pn] for pn in all_pn)
        init_states = gen_init_states(all_p)
        
        for i in 1:length(all_vn)
            prob.u0[i] = init_states[all_vn[i]]
        end
            
        return remake(prob,u0=convert.(eltype(p),prob.u0),p=p)
    end
end # get_prob_update

"""
get_cost_function(fitted_vn::Vector{String},    # variable names that are fitted
    all_vn::Vector{String},       # variable names: gives the order of variables
    time_points::Vector{Float64}, # time points of data values to fit
    prob_update::Function,        # problem updater of init states
    data::Array{Float64,2},       # time points of data values to fit
    prob::ODEProblem)            # problem

returns the cost function for single set of data points.
"""
function get_cost_function(fitted_vn::Vector{String},    # variable names that are fitted
                           all_vn::Vector{String},       # variable names: gives the order of variables
                           time_points::Vector{Float64}, # time points of data values to fit
                           prob_update::Function,        # problem updater of init states
                           data::Array{Float64,2},       # time points of data values to fit
                           prob::ODEProblem)            # problem
    build_loss_objective(prob,Tsit5(),L2Loss(time_points,data),prob_generator = prob_update, maxiters=1000000,verbose = false, save_idxs = indexin(fitted_vn,all_vn))
end

"""
get_cost_function(fitted_vn::Vector{String},    # variable names that are fitted
    all_vn::Vector{String},       # variable names: gives the order of variables
    time_points::Vector{Float64}, # time points of data values to fit
    prob_update::Function,        # problem updater of init states
    data::Matrix{Normal{Float64}},# distribution of data values to fit
    prob::ODEProblem)            # problem

returns the cost function for distributions of data points.
"""
function get_cost_function(fitted_vn::Vector{String},    # variable names that are fitted
                           all_vn::Vector{String},       # variable names: gives the order of variables
                           time_points::Vector{Float64}, # time points of data values to fit
                           prob_update::Function,        # problem updater of init states
                           data::Matrix{Normal{Float64}},# distribution of data values to fit
                           prob::ODEProblem)            # problem
    build_loss_objective(prob,Tsit5(),LogLikeLoss(time_points,data),prob_generator = prob_update, maxiters=1000000,verbose = false, save_idxs = indexin(fitted_vn,all_vn))
end

"""
param_estim(
    free_pn::Vector{String},      # free parameter names
    all_p::Dict{String,Float64},  # parameter -> value dict
    fitted_vn::Vector{String},    # variable names that are fitted
    all_vn::Vector{String},       # variable names: gives the order of variables
    lower::Dict{String,Float64},  # lower bounds for parameter search
    upper::Dict{String,Float64},  # upper bounds for parameter search
    time_points::Vector{Float64}, # time points of data values to fit
    data,                         # data values to fit
    prob::ODEProblem)

performs the parameter estimation for the ODE problem `prob` of parameters listed in `free_pn`.
"""
function param_estim(
                     free_pn::Vector{String},      # free parameter names
                     all_p::Dict{String,Float64},  # parameter -> value dict
                     fitted_vn::Vector{String},    # variable names that are fitted
                     all_vn::Vector{String},       # variable names: gives the order of variables
                     lower::Dict{String,Float64},  # lower bounds for parameter search
                     upper::Dict{String,Float64},  # upper bounds for parameter search
                     time_points::Vector{Float64}, # time points of data values to fit
                     data,                         # data values to fit
                     prob::ODEProblem,             # problem
                     my_optimize                   # optimization function
                     )

    free_p = [ all_p[pn] for pn in free_pn ] # vector of seeds for free parameter values
    
    free_lower::Vector{Float64} = [ lower[pn] for pn in free_pn ]
    free_upper::Vector{Float64} = [ upper[pn] for pn in free_pn ]

    prob_update = get_prob_update(free_pn, all_p, all_vn)
    cost_function = get_cost_function(fitted_vn,all_vn,time_points,prob_update,data,prob)
    return my_optimize(cost_function, free_lower, free_upper, free_p)
end # param_estim

"""
    update_params(free_pn::Vector{String}, 
                  free_pv::Vector{Float64}, 
                  all_pn::Vector{String}, 
                  all_p::Dict{String,Float64})::Dict{String,Float64}
    returns a dictionary of parameters values of all_p updated with the values in free_pv for the parameters of free_pn.
"""
function update_params(free_pn::Vector{String}, 
                       free_pv::Vector{Float64}, 
                       all_pn::Vector{String}, 
                       all_p::Dict{String,Float64})::Dict{String,Float64}
    idx_of_free_pn = Dict(pn => i for (pn,i) in zip(free_pn,1:length(free_pn)))
    Dict(pn => (pn in free_pn) ? free_pv[idx_of_free_pn[pn]] : all_p[pn] for pn in all_pn)
end # update_params
