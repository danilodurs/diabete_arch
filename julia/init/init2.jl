description = "fit of Gp without basal values and merged compartments"

# variables used for curve fitting
fitted_vn = ["Gp"]

# time span interval of integration
tspan = (0.0,420.0)

# compartments for parameter estimation with basal estimation
compartments = Dict(
     "GK+U" => ["pk_1", "pk_2", "pV_m0", "pV_mx", "pxK_m0", "pp_2U"],
     "IK+S" => ["pV_I", "pm_1", "pm_5", "pm_6", "pK", "pxα_Y", "pxβ_Y"], # we assume pHEb constant
     "Ra"   => ["pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc"],
     "EGP"   => ["pk_p2", "pk_p3", "pk_p4", "pk_i"],
     "RE"   => ["pk_e1", "pk_e2"],
)

# coefficient used to adjust the width of the parameter interval search
coeff_bounds = 2.0

# delay between to simulated data points
delay = 30.0
