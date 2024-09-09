using BlackBoxOptim

function my_optimize(cost_function, free_lower, free_upper, free_p)
    #BlackBoxOptim.bboptimize(cost_function,free_p;SearchRange = collect(zip(free_lower,free_upper)), MaxSteps = 11e3)
    BlackBoxOptim.bboptimize(cost_function;SearchRange = collect(zip(free_lower,free_upper)), MaxSteps = 11e3)
end

function best_params(res)
    best_candidate(res)
end

function loss_value(res)
    best_fitness(res)
end