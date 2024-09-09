using Optim

function my_optimize(cost_function, free_lower, free_upper, free_p)
    Optim.optimize(cost_function, free_lower, free_upper, free_p, Fminbox(BFGS()))
end 

function best_params(res)
    res.minimizer
end

function loss_value(res)
    res.minimum
end