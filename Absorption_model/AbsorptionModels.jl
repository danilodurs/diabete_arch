module AbsorptionModels

using DifferentialEquations
using Catalyst
using Plots
using DiffEqParamEstim
using Optim

export gen_model, gen_prob_partiel, gaussian, discretiz
export gen_data_set, gen_dist_transit, gen_Ra

"""
        gen_model(n)

Generate a reaction network model of absorption with `n` compartments.
"""
function gen_model(n)
        if (n<1)
                throw(DomainError(n,"Cannot generate model with less than 1 compartment!!"))
        end

        # pure mass action (and thus exponential) kinetics of gastric emptying
        # Qgutn = Symbol("Qgut", n)
        # rn = @eval @reaction_network begin
        #        kempt, Qsto → Qgut1
        #        kabs1, Qgut1 → Qp
        #        kex, $Qgutn → ∅
        # end kempt kabs1 kex

        # gastric emptying as defined in Salinari et al.
        @reaction_func gastric_empt(D,k) = D*1.2*k^(1.2)*t^(.2)*exp(-(k*t)^(1.2))

        Qgutn = Symbol("Qgut", n)
        rn = @eval @reaction_network begin
                gastric_empt(416.3,kempt), Qsto ⇒ Qgut1
                kabs1, Qgut1 → Qp
                kex, $Qgutn → ∅
        end kempt kabs1 kex

        for i in 1:n-1
                Qgut = Symbol("Qgut",i)
                QgutNext = Symbol("Qgut",i+1)
                kk = Symbol("k",i,i+1)
                kabs = Symbol("kabs",i+1)
                @eval merge!($rn,@reaction_network begin
                       $kk, $Qgut → $QgutNext
                       $kabs, $QgutNext → Qp
               end $kk $kabs)
       end
       rn
end

"""
        gen_prob_partiel(model, param_values, tspan, u0)

Compute a ODE problem from a reaction network `model`, a (possibly partial)
dictionary `param_values` of parameter/values, time span `tspan` and initial
condition `u0`.
"""
function gen_prob_partiel(model, param_values, tspan, u0)
        ode_fun = ODEFunction(convert(ODESystem,model))
        the_parameters = params(model)
        nb_params = numparams(model)
        function ode_fun_part!(du,u,p,t)
                pp = zeros(nb_params)
                j = 1
                for i in 1:nb_params
                        if the_parameters[i] in keys(param_values)
                                pp[i] = param_values[the_parameters[i]]
                        else
                                pp[i] = p[j]
                                j += 1
                        end
                end
                ode_fun(du,u,pp,t)
        end
        ODEProblem((du,u,p,t) -> Base.invokelatest(ode_fun_part!,du,u,p,t),u0,tspan,param_values)
end

"""
        gaussian(μ,σ,c)

Compute of a Gaussian function with mean μ, standard deviation σ and maximum
inversely proportional ro `c`.
"""
function gaussian(μ,σ,c)
        return x -> (c/(σ*sqrt(2*pi)))*exp(-0.5*((x-μ)^2/σ^2))
end

"""
        discretiz(f,a,b,n)

Output a discretization of `f` in the interval [`a`,`b`] with `n` intervals.
"""
function discretiz(f,a,b,n)
        δ = (b-a)/n # width of one compartment
        p = δ # precision of sampling
        result = Vector()
        for i in 0:n-1
                range = i*δ:δ/p:(i+(p-1)/p)*δ
                push!(result,sum(f.(range))/p)
        end
        result
end

"""
        gen_data_set(model,param_values,tspan,u0)

Generates a data set for a given model, parameter values, integration interval
and initial conditions.
"""
function gen_data_set(model,param_values,tspan,u0)
       prob = gen_prob_partiel(model, param_values, tspan, u0)
       solve(prob, Tsit5())
end

"""
        gen_data_set(n, abs_dist, trans_dist, kempt_v, kex_v)

Generates a data set for a given number of compartments `n`, distribution
function of absorption `abs_dist`, and distribution function of intestinal
transit rates `trans_rate`, gastric emptying rate `kempt_v` and excretion
rate `kex_v`.
"""
function gen_data_set(n, abs_dist, trans_dist, kempt_v, kex_v)
        # generate a model of n compartments
        model = gen_model(n)

        # infer parameter identifiers
        kempt, kabs1, kex = params(model)[1:3]
        for i in 1:n-1
               local k12 = Symbol("k",i,i+1)
               local kabs2 = Symbol("kabs",i+1)
               @eval $k12,$kabs2 = params($model)[2*$i+2:2*$i+3]
        end

        # set the initial state u0 and the integration interval tspan
        # The initial glucose amount is 75g that is 416.3 mmol
        u0 = [416.3, .0, .0]
        for i in 2:n
                push!(u0,.0)
        end
        tspan = (0.0, 600.0)

        # set the parameter values
        kabs_v = discretiz(abs_dist, 0, 600, n)
        ktrans_v = discretiz(trans_dist, 0, 600, n-1)
        param_values = Dict(kempt => kempt_v, kabs1 => kabs_v[1], kex => kex_v)
        for i in 1:n-1
                local k12 = Symbol("k", i, i+1)
                local kabs2 = Symbol("kabs", i+1)
                @eval $param_values[$k12] = $ktrans_v[$i]
                @eval $param_values[$kabs2] = $kabs_v[$i+1]
        end

        # data simulation
        gen_data_set(model,param_values,tspan,u0)
end

"""
        gen_dist_transit(sol, t, n)

Given the simulation `sol` an absorption model, return the distribution of the glucose concentration
in the intestin at time `t` as a vector of length `n`.
"""
function gen_dist_transit(sol, t, n)
    dist_trans = [sol(t)[2]]
    for i in 5:n+2
            push!(dist_trans, sol(t)[i])
    end
    if (n>1)
            push!(dist_trans, sol(t)[4])
    end
    dist_trans
end

"""
        gen_Ra(sol, abs_dist, n)

Plot the rate of absorption for given data set solution `sol`, absorption
distribution `abs_dist` and number of compartments n
"""
function gen_Ra(sol, abs_dist, n)
        kabs_v = discretiz(abs_dist, 0, 600, n)
        Ra = kabs_v[1]*sol[2,:]
        for i in 5:n+2
                Ra = Ra + kabs_v[i-3]*sol[i,:]
        end
        if n > 1
                Ra = Ra + kabs_v[n]*sol[4,:]
        end
end

end  # module AbsorptionModels
