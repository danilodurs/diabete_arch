using CMAEvolutionStrategy

function my_optimize(cost_function, free_lower, free_upper, free_p)
    minimize(cost_function,free_p,1.;lower=free_lower,upper=free_upper)
end

function best_params(res)
    xbest(res)
end

function loss_value(res)
    fbest(res)
end