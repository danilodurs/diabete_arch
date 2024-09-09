

##################################
### lower_upper_bounds function  #
##################################

# TODO: modifier cette fonction pour qu'elle renvoie des intervalles en échelle linéaire ou logarithmique, c'est-à-dire:
# étant données 2 valeurs de paramètres p1 et p2 (typiquement, les valeurs de patient normal et diabétique) telles que 
# p2 = k*p1 (p2 est k fois plus grand que p1), alors:
# - soit k ≥ 100 et alors on fixe la borne inf p_min = p1/k et la borne sup p_max = p2*k,
# - soit k < 100 et alors on fixe la borne inf p_min = max(0,p1-Δp) et la borne sup p_max = p2+Δp
# où Δp = coef*|p1-p2| et où coef vaut, par défaut, 10.
function lower_upper_bounds(params::Dict{Symbol,Vector{Float64}}, # maps each parameter name to pairs of parameters values
    coef=10::Float64;                   # determines the width of upper and lower bound
    st=1                           # specifies around which subject type parameters the interval have to be built (1 or 2)
)::Tuple{Dict{Symbol,Float64},Dict{Symbol,Float64}}
    params_k = Dict(pn => (max(params[pn][2], params[pn][1]) / min(params[pn][2], params[pn][1])) for pn in keys(params))
    params_Δ = Dict(pn => coef * abs(params[pn][2] - params[pn][1]) for pn in keys(params))
    lower = Dict{Symbol,Float64}()
    upper = Dict{Symbol,Float64}()
    for pn in keys(params)
        if (params_k[pn] < 100)
            lower[pn] = max(0, min(params[pn][1], params[pn][2]) - params_Δ[pn])
            upper[pn] = max(params[pn][1], params[pn][2]) + params_Δ[pn]
        else
            lower[pn] = min(params[pn][1], params[pn][2]) / params_k[pn]
            upper[pn] = max(params[pn][1], params[pn][2]) * params_k[pn]
        end
    end
    (lower, upper)
end

# CAUTION: this function does not check if some parameters have a lower bound equal to an upper bound
"""
    lower_upper_bounds(params::Dict{String,Vector{Float64}}, # maps each parameter name to pairs of parameters values
                       coef::Float64;                        # determines the width of upper and lower bound
                       st = "diab"                           # specifies around which subject type parameters the interval have to be built
                      )::Tuple{Vector{Float64},Vector{Float64}}
"""
function lower_upper_bounds(params::Dict{Symbol,Vector{Float64}}, # maps each parameter name to pairs of parameters values
    coef::Float64;                        # determines the width of upper and lower bound
    st=1                           # specifies around which subject type parameters the interval have to be built (1 or 2)
)::Tuple{Dict{Symbol,Float64},Dict{Symbol,Float64}}
    params_Δ = Dict(pn => coef * abs(params[pn][2] - params[pn][1]) for pn in keys(params))
    lower = Dict(pn => max(0, params[pn][st] - params_Δ[pn]) for pn in keys(params))
    upper = Dict(pn => params[pn][st] + params_Δ[pn] for pn in keys(params))
    # if the other subject type's parameters are not reachable we raise an exception
    for pn in keys(params)
        if !(lower[pn] <= params[pn][(st%2)+1] <= upper[pn])
            throw(DomainError(pn, "normal parameters are not reachable!"))
        end
    end
    (lower, upper)
end # lower_upper_bounds

##########################
### Parameter estimation #
##########################

# usefull functions

# My cost function (without normalization which is useless because we consider only one data time series) 
# Useless function 
function gen_cost_function(prob, timeticks, data)
    function my_cost_function(prm)
        prob = remake(prob; p=prm)
        sol = solve(prob)
        loss = 0
        for (t, i) in zip(timeticks, 1:length(timeticks))
            loss += (sol(t)[1] - data[i])^2
        end
        loss
    end
    my_cost_function
end

# functions for parameter normalization (this assumes that xmin and xmax != 0 other it always returns 0!) 
function from_log_010_to_ab(x, xmax, xmin)
    return xmin * (xmax / xmin)^(x / 10)
end

"""
    prepare_estimation(crn, u0, timeticks, data)

Returns cost function, initial parameter values and upper/lower parameter bounds 
"""
function prepare_estimation(crn, u0, timeticks, data)
    ode_model = convert(ODESystem, crn)

    p0 = [all_prm[prm_symbol][1] for prm_symbol in Symbol.(parameters(crn))]
    (all_lower, all_upper) = lower_upper_bounds(all_prm, 10.0)
    lowerb = [all_lower[prm_symbol] for prm_symbol in Symbol.(parameters(crn))]
    upperb = [all_upper[prm_symbol] for prm_symbol in Symbol.(parameters(crn))]

    tspan = (0.0, timeticks[end])
    prob = ODEProblem(ode_model, u0, tspan, p0)

    cost_function = build_loss_objective(prob, Tsit5(), L2Loss(timeticks, data), maxiters=1000000, verbose=false, save_idxs=1)

    (cost_function, p0, lowerb, upperb)
end

# for given parameter values p0, cost_function(p0) returns the loss value.
# We here provide a wrapper wrapped_cost_function that will take normalized 
# parameter values as arguments and call cost_function on its "real" parameter 
# values
function build_wrapped_cost_function(crn, u0, timeticks, data)
    function my_cost_function(p0)
        loss = cost_function(from_log_010_ab())
    end
end

# run just many times the param estimation for a given crn model. Outputs min parameters and loss.
function run_estimation(cost_function, p0, lowerb, upperb, max_run)
    min_prm = p0
    min_loss = 10^10
    for i in range(1, max_run)
        res = my_optimize(cost_function, lowerb, upperb, p0)
        # res = minimize(cost_function,p0,1.,lower = lowerb, upper = upperb)
        if (loss_value(res) < min_loss)
            min_loss = loss_value(res)
            min_prm = best_params(res)
        end
    end
    (min_prm, min_loss)
end

"""
    run_all_estimations(crn, u0, idstr, timeticks, data)

Run parameter estimation for each minipig
"""
function run_all_estimations(crn, u0, idstr, timeticks, data)
    results = Dict()
    for minipig_id in idstr
        data_vec = Matrix{Float64}(data[minipig_id]')
        u0[1] = data_vec[1]
        (cost_function, p0, lowerb, upperb) = prepare_estimation(crn, u0, timeticks, data_vec)
        results[minipig_id] = run_estimation(cost_function, p0, lowerb, upperb, 10)
    end
    results
end

"""
    run_all_estimations(idstr, timeticks, data)

version specific to power exponential model
"""
function run_all_estimations(idstr, timeticks, data)
    results = Dict()
    for minipig_id in idstr
        data_vec = Matrix{Float64}(data[minipig_id]')
        DXp0 = maximum(data_vec)
        (cost_function, p0, lowerb, upperb) = prepare_estimation(gen_crn_powexp(DXp0), [DXp0], timeticks, data_vec)
        results[minipig_id] = run_estimation(cost_function, p0, lowerb, upperb, 10)
    end
    results
end

# run a simulation for a given ODEProblem, returns the plot displayed with the given fig name, loss and parameter values 
function plot_simulation(prob, prm, u0, data, timeticks, loss, fig_name)
    u0[1] = data[1]
    prob = remake(prob; p=prm, u0=u0)
    sol = solve(prob)
    plt = plot(sol, xlabel="t", ylabel="mg/kg", label="DXp", vars=(0, 1), title=fig_name * " (loss = " * string(round(loss)) * ")")
    plt = plot!(timeticks, data, m=:circle, label="DX data")
    plt = annotate!(580, trunc(Int, maximum(data)) - 100, text("prm = " * string(round.(prm, sigdigits=3)), :grey, :right, 10))
    plt
end

# make plots of all parameter estimation results (given as a dictionary), returns all the plots as dictionary with minipig IDs for keys.
function make_plots(crn, u0, results, data, timeticks)
    tspan = (0.0, timeticks[end])
    prob = ODEProblem(convert(ODESystem, crn), u0, tspan, zeros(length(parameters(crn))))
    plts = Dict()
    for minipig_id in keys(results)
        plts[minipig_id] = plot_simulation(prob, results[minipig_id][1], u0, data[minipig_id], timeticks, results[minipig_id][2], minipig_id)
    end
    plts
end

# specific version for power exponential model
function make_plots(results, data, timeticks)
    tspan = (0.0, timeticks[end])
    plts = Dict()
    for minipig_id in keys(results)
        DXp0 = maximum(data[minipig_id])
        crn = gen_crn_powexp(DXp0)
        prob = ODEProblem(convert(ODESystem, crn), [DXp0], tspan, zeros(length(parameters(crn))))
        plts[minipig_id] = plot_simulation(prob, results[minipig_id][1], [DXp0], data[minipig_id], timeticks, results[minipig_id][2], minipig_id)
    end
    plts
end

# save all the plots (given as dictionary with minipig IDs for keys) in a given directory.
function save_plots(plts, directory_name)
    for minipig_id in keys(plts)
        savefig(plts[minipig_id], directory_name * "/" * minipig_id * ".png")
    end
end

