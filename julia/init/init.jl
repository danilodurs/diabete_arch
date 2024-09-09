module Init

export fitted_vn
export tspan
export compartments
export coeff_bounds
export delay

# fit of Gp and basal values and merged compartments
# include("./init1.jl")
# fit of Gp without basal values and merged compartments
include("./init2.jl")

end # Init
