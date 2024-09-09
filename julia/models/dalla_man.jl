
# module Dalla_Man

# export Sb, m_4, m_2, Gtb, V_m0, k_p1
# export dalla_net!
# export params
# export all_pn, all_p
# export all_vn, init_v


# computation of some parameters
function Sb(pm_5,pm_6,pHE_b) # from eq (6)
    (pm_6-pHE_b)/pm_5
end

function m_4(pm_5,pm_6,pIb,pV_I,pHE_b) # from eq (9)
    (2/5)*(Sb(pm_5,pm_6,pHE_b)/(pIb*pV_I))*(1-pHE_b)
end

function m_2(pm_5,pm_6,pIb,pV_I,pHE_b) # from eq (9)
    (Sb(pm_5,pm_6,pHE_b)/(pIb*pV_I)-(m_4(pm_5,pm_6,pIb,pV_I,pHE_b)/(1-pHE_b)))*(1-pHE_b)/pHE_b
end

# to compute the basal glucose Gb, one should solve the following equations:
# using Reduce
# eq1 = :(pEGPb = Uidb + pF_cns)
# eq2 = :(Gtb * pV_m0 = (pEGPb-pF_cns)*(pxK_m0 + Gtb))
# eq3 = :(pk_1 * Gpb = pEGPb - pF_cns + pk_2*Gtb)
# sol = Algebra.solve((eq1,eq2,eq3),(:Gpb,:Uidb,:Gtb))
# sol[1]

function Gb(pk_1, pk_2, pV_m0, pxK_m0, pF_cns, pEGPb, pV_G)
    ((pk_2 * pxK_m0 + pV_m0 + 2pF_cns) * pEGPb - ((pk_2 * pxK_m0 + pV_m0 + pF_cns) * pF_cns + pEGPb ^ 2)) / (((pF_cns + pV_m0) - pEGPb) * pk_1)/pV_G
end

# function Gtb(pF_cns,pEGPb,pk_1,pV_G,pk_2,pV_m0,pxK_m0) # from eq (20)
#     pGb = Gb(pk_1, pk_2, pV_m0, pxK_m0, pF_cns, pEGPb)
#     (pF_cns-pEGPb+pk_1*pGb*pV_G)/pk_2
# end

# function V_m0(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2,pxK_m0) # from eq (22)
#     (pEGPb-pF_cns)*(pxK_m0+Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2))/Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2)
# end

function k_p1(pk_p2,pk_p3,pk_p4,pEGPb,pGb,pV_G,pγ,pSb,pIb) # from eq (12)
    pEGPb + pk_p2*pGb*pV_G+pk_p3*pIb+pk_p4*pSb/pγ
end

# model definition

function dalla_net!(du, u, p, t)
    # parameters
    pV_G, pk_1, pk_2, pV_I, pm_1, pm_5, pm_6, pHE_b, pk_max, pk_min, pk_abs, pk_gri, pf, pb, pc, pBW, pD, pDX, pk_p2, pk_p3, pk_p4, pk_i, pF_cns, pV_m0, pV_mx, pxK_m0, pxK_mx, pp_2U, pK, pxα_Y, pxβ_Y, pγ, pk_e1, pk_e2, DXpk_e1, DXpk_e2, pIb, pEGPb = p

    pSb = Sb(pm_5,pm_6,pHE_b)
    # pGtb = Gtb(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2)
    # pV_m0 = V_m0(pF_cns,pEGPb,pk_1,pGb,pV_G,pk_2,pxK_m0)
    pGb = Gb(pk_1, pk_2, pV_m0, pxK_m0, pF_cns, pEGPb, pV_G)
    pk_p1 = k_p1(pk_p2,pk_p3,pk_p4,pEGPb,pGb,pV_G,pγ,pSb,pIb)
    pm_4 = m_4(pm_5,pm_6,pIb,pV_I,pHE_b)
    pm_2 = m_2(pm_5,pm_6,pIb,pV_I,pHE_b)
    pα = 5 / (2*pD*(1-pb))
    pβ = 5 / (2*pc*pD)
    ph = pGb

    # algebraic equations
    G     = t -> u[1] / pV_G # G

    I     = t -> u[4] / pV_I                  # I
    m_3   = t -> (HE(t) * pm_1) / (1 - HE(t)) # m_3
    HE    = t -> - pm_5 * S(t) + pm_6          # HE

    Ra     = t     -> (pf * pk_abs * u[7]) / pBW # Ra
    q_sto  = t     -> u[5] + u[6]             # q_sto
    k_empt = q_sto -> pk_min + (( pk_max -  pk_min) / 2.) * (tanh( pα * (q_sto -  pb *  pD)) - tanh( pβ * (q_sto - pc *  pD))+ 2.) # k_empt 

    EGP   = t -> max(pk_p1 - pk_p2 * u[1] - pk_p3 * u[9] - pk_p4 * u[12],0) # EGP

    U_ii   = t   -> pF_cns                                     # U_ii
    U_id   = t   -> (V_m(u[10]) * u[2]) / (K_m(u[10]) + u[2]) # U_id
    V_m    = X_t -> pV_m0 + pV_mx * X_t                         # V_m
    K_m    = X_t -> pxK_m0 + pxK_mx * X_t                         # K_m

    function piecewiseDY(t)          # d.Y
        if pxβ_Y * (G(t) - ph) >= - pSb
            return - pxα_Y * (u[11] - pxβ_Y * (G(t) - ph))
        elseif pxβ_Y * (G(t) - ph) < - pSb
            return - pxα_Y * u[11] - pxα_Y * pSb
        end
    end
    S      = t -> pγ * u[12]          # S
    S_po   = t ->  piecewiseS_po(t)  # S_po
    function piecewiseS_po(t)        # S_po   
        if (du[1] / pV_G) > 0.
            return u[11] + pK * (du[1] / pV_G) + pSb
        elseif (du[1]/ pV_G) <= 0.
            return u[11] + pSb
        end
    end
    # ---------------------
    # Renal Excretion
    E = t -> piecewiseE(t) # E
    function piecewiseE(t) # E
        if u[1] > pk_e2
            return pk_e1 * (u[1] - pk_e2)
        elseif u[1] <= pk_e2
            return 0.
        end
    end
    
    # D-Xylose algebraic equations
    
    DXRa     = t     -> (pf * pk_abs * u[17]) / pBW # DXRa
    DXq_sto  = t     -> u[15] + u[16]             # DXq_sto
    DXk_empt = DXq_sto -> pk_min + (( pk_max -  pk_min) / 2.0) * (tanh( pα * (DXq_sto -  pb *  pDX)) - tanh( pβ * (DXq_sto - pc *  pDX))+ 2.0) # DXk_empt 
    # D-Xylose Renal Excretion
    DXE = t -> piecewiseDXE(t) # DXE
    function piecewiseDXE(t) # E
        if u[13] > DXpk_e2
            return DXpk_e1 * (u[13] - DXpk_e2)
        elseif u[13] <= DXpk_e2
            return 0.0
        end
    end
    
    # ------------------------------------------------------------------------
    # differential equations

    # Glucose kinetics
    du[1] = EGP(t) + Ra(t) - U_ii(t) - E(t) - pk_1 * u[1] + pk_2 * u[2] # d.G_p
    du[2] = - U_id(t) + pk_1 * u[1] - pk_2 * u[2]                       # d.G_t


    # --------------------
    # Insulin kinetics
    du[3] = - (pm_1 + m_3(t)) * u[3] + pm_2 * u[4] + S(t) # d.I_l
    du[4] = - (pm_2 + pm_4) * u[4] + pm_1 * u[3]           # d.I_p


    # -------------------- 
    # Glucose rate of appearance
    #du[5]  = - pk_gri * u[5] + pD * DiracDelta(0.001)(t) # d.q_sto1
    du[5]  = - pk_gri * u[5] # d.q_sto1
    du[6]  = - k_empt(q_sto(t)) * u[6] + pk_gri * u[5]  # d.q_sto2
    du[7]  = - pk_abs * u[7] + k_empt(q_sto(t)) * u[6]  # d.q_gut


    # --------------------
    # Endogenous production
    du[8] = - pk_i * (u[8] - I(t)) # d.I_1
    du[9] = - pk_i * (u[9] - u[8]) # d.I_d


    # --------------------
    # Utilization
    du[10] = - pp_2U * u[10] + pp_2U * (I(t)-pIb) # d.X


    # --------------------
    # Secretion
    du[11] = piecewiseDY(t)          # d.Y
    du[12] = - pγ * u[12] + S_po(t)   # d.I_po   
    
    # D-Xylose ODEs-------------------------------------
    # D-Xylose kinetics
     du[13] =  DXRa(t) - DXE(t) - pk_1 * u[13] + pk_2 * u[14] # d.DX_p
    du[14] = pk_1 * u[13] - pk_2 * u[14]                       # d.DX_t
    # D-Xylose Gastro-Intestinal System
    du[15] = - pk_gri * u[15] # d.q_sto1
    du[16] = - DXk_empt(DXq_sto(t)) * u[16] + pk_gri * u[15]  # d.DXq_sto2
    du[17] = - pk_abs * u[17] + DXk_empt(DXq_sto(t)) * u[16]  # d.DXq_gut
    
    
end # dalla_net!

# parameter values (both normal and diabetic)

DM_params = Dict("pV_G"   => [1.88,1.49], 
              "pk_1"   => [0.065,0.042],  
              "pk_2"   => [0.079,0.071], 
              "pV_I"   => [0.05,0.04], 
              "pm_1"   => [0.190,0.379], 
              # "pm_2"   => [0.484,0.637], # can be computed from equation (9)
              # "pm_4"   => [0.194,0.269], # can be computed from equation (9)
              "pm_5"   => [0.0304,0.0526], 
              "pm_6"   => [0.6471,0.8118],
              "pHE_b"  => [0.6,0.6],
              "pk_max" => [0.0558,0.0465], 
              "pk_min" => [0.0080,0.0076], 
              "pk_abs" => [0.057,0.023], 
              "pk_gri" => [0.0558,0.0465],  
              "pf"     => [0.9,0.9], 
              "pb"     => [0.82,0.68], 
              "pc"     => [0.010,0.09], 
              # "pα"     => [0.00013,0.00006], # can be computed from eq (10) in [DallaMan & al. 2006]
              # "pβ"     => [0.00236,0.00023], # can be computed from eq (11) in [DallaMan & al. 2006]
              "pBW"    => [78.0,91.0],
              "pD"     => [90000.0, 90000.0], 
              "pDX"    => [30000.0,30000.0],
              # "pk_p1"  => [2.70,3.09], # can be computed from equation (12)
              "pk_p2"  => [0.0021,0.0007], 
              "pk_p3"  => [0.009,0.005], 
              "pk_p4"  => [0.0618,0.0786], 
              "pk_i"   => [0.0079,0.0066], 
              "pF_cns" => [1.0,1.0],
              "pV_m0"  => [2.50,4.65], # We keep this parameter even if it can be computed as in eq 22 of Dalla Man's paper
              "pV_mx"  => [0.047,0.034], 
              "pxK_m0" => [225.59,466.21], 
              "pxK_mx" => [0.0,0.0],
              "pp_2U"  => [0.0331,0.0840], 
              "pK"     => [2.30,0.99], 
              "pxα_Y"  => [0.050,0.013], 
              "pxβ_Y"  => [0.11,0.05], 
              "pγ"     => [0.5,0.5], 
              # "ph"     => [89.39980053,141.32659636], # h is supposed = to Gb
              "pk_e1"  => [0.0005,0.0007], 
              "pk_e2"  => [339.0,269.0],
              "DXpk_e1"  => [0.0005,0.0007], 
              "DXpk_e2"  => [0.001, 0.001],
              "pIb"    => [25.49, 54.81], 
              # "pSb"    => [1.54,3.57], # I remove this parameter because it is computed from HEb: pSb = (m6-HEb)/m5
              # we add 2 parameters to allow for estimation of basal Glucose and EGP
              # "pGb"    => [91.76, 164.18],
              "pEGPb"  => [1.92, 2.01])

all_pn = ["pV_G", "pk_1", "pk_2", "pV_I", "pm_1", "pm_5", "pm_6", "pHE_b", "pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc", "pBW", "pD", "pDX", "pk_p2", "pk_p3", "pk_p4", "pk_i", "pF_cns", "pV_m0", "pV_mx", "pxK_m0", "pxK_mx", "pp_2U", "pK", "pxα_Y", "pxβ_Y", "pγ", "pk_e1", "pk_e2", "DXpk_e1", "DXpk_e2", "pIb", "pEGPb"]
#all_pn = ["pV_G", "pk_1", "pk_2", "pV_I", "pm_1", "pm_5", "pm_6", "pHE_b", "pk_max", "pk_min", "pk_abs", "pk_gri", "pf", "pb", "pc", "pBW", "pD", "pDX", "pk_p2", "pk_p3", "pk_p4", "pk_i", "pF_cns", "pV_m0", "pV_mx", "pxK_m0", "pxK_mx", "pp_2U", "pK", "pxα_Y", "pxβ_Y", "pγ", "pk_e1", "pk_e2", "pIb", "pEGPb"]
# parameters are more convenient to handle with the following dictionary
all_p = Dict(st => Dict(pn => DM_params[pn][i] for pn in all_pn) 
             for (st,i) in zip(["norm","diab"],[1,2]))

# all variable names
all_vn = ["Gp","Gt","Il","Ip","Qsto1","Qsto2","Qgut","I1","Id","X","Y","Ipo","DXp","DXt","DXQsto1","DXQsto2","DXQgut"]

# end # Dalla_Man module

"""
    traces_of_obs_vars(sol,                         # solution of an ODE
                       all_p::Dict{String,Float64}  # all parameter values
                      )::Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}

return the tuple of the traces of G, I, EGP, Ra, U and S.
"""

function traces_of_obs_vars(sol,all_p)#::Dict{String,Float64})::Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}
    
    # pGtb = Gtb(all_p["pF_cns"],all_p["pEGPb"],all_p["pk_1"],all_p["pGb"],all_p["pV_G"],all_p["pk_2"])
    #pV_m0 = V_m0(all_p["pF_cns"],all_p["pEGPb"],all_p["pk_1"],all_p["pGb"],all_p["pV_G"],all_p["pk_2"],all_p["pxK_m0"])
    pSb = Sb(all_p["pm_5"],all_p["pm_6"],all_p["pHE_b"])
    pGb = Gb(all_p["pk_1"],all_p["pk_2"], all_p["pV_m0"], all_p["pxK_m0"], all_p["pF_cns"], all_p["pEGPb"], all_p["pV_G"])
    pk_p1 = k_p1(all_p["pk_p2"],all_p["pk_p3"],all_p["pk_p4"],all_p["pEGPb"],pGb,all_p["pV_G"],all_p["pγ"],pSb,all_p["pIb"])

    Vm = X -> all_p["pV_m0"] + all_p["pV_mx"] * X
    Km = X -> all_p["pxK_m0"] + all_p["pxK_mx"] * X

    G = sol[1,:] / all_p["pV_G"] 
    I = sol[4,:] / all_p["pV_I"] 
    EGP = -all_p["pk_p2"]*sol[1,:] - all_p["pk_p3"]*sol[9,:] - all_p["pk_p4"]*sol[12,:] .+ pk_p1
    Ra = all_p["pf"] * all_p["pk_abs"] * sol[7,:] / all_p["pBW"]

    DXp = sol[13,:]
    U = (Vm.(sol[10,:]) .* sol[2,:])./(Km.(sol[10,:]) .+ sol[2,:]) .+ all_p["pF_cns"]
    S = all_p["pγ"] * sol[12,:]

    (G,I,EGP,Ra,DXp,U,S)


end # traces_of_obs_vars

"""
    gen_init_states(params::Dict{String,Float64})::Dict{String,Float64}

returns the init states of all variables as a dictionary
"""
function gen_init_states(params::Dict{String,Float64})::Dict{String,Float64}
    init_states::Dict{String,Float64} = Dict()
    m30::Float64 = params["pHE_b"]*params["pm_1"]/(1-params["pHE_b"])
    pSb = Sb(params["pm_5"],params["pm_6"],params["pHE_b"])
    pm_4 = m_4(params["pm_5"],params["pm_6"],params["pIb"],params["pV_I"],params["pHE_b"])
    pm_2 = m_2(params["pm_5"],params["pm_6"],params["pIb"],params["pV_I"],params["pHE_b"])
    pGb = Gb(params["pk_1"],params["pk_2"], params["pV_m0"], params["pxK_m0"], params["pF_cns"], params["pEGPb"], params["pV_G"])

    # pSb = (params["pm_6"]-params["pHE_b"])/params["pm_5"]
    # pm_4 = (2/5)*(pSb/(params["pIb"]*params["pV_I"]))*(1-params["pHE_b"]) # from eq (9)
    # pm_2 = (pSb/(params["pIb"]*params["pV_I"])-(pm_4/(1-params["pHE_b"])))*(1-params["pHE_b"])/params["pHE_b"] # from eq (9)

    init_states["Gp"] = pGb*params["pV_G"]
    init_states["Gt"] = (params["pF_cns"]-params["pEGPb"]+params["pk_1"]*init_states["Gp"])/params["pk_2"]
    init_states["Ip"] = params["pIb"]*params["pV_I"]
    init_states["Il"] = init_states["Ip"]*(pm_4+pm_2)/params["pm_1"]
    init_states["Qsto1"] = params["pD"]
    init_states["Qsto2"] = 0.0
    init_states["Qgut"] = 0.0
    init_states["I1"] = params["pIb"]
    init_states["Id"] = params["pIb"]
    init_states["X"] = 0.0
    init_states["Y"] = 0.0
    #init_states["Ipo"] = (init_states["Ip"]*pm_4+init_states["Il"]*m30)/params["pγ"]
    init_states["Ipo"] = pSb/params["pγ"] # see eq (23) in Dalla Man's paper

    init_states["DXp"] = 0.0
    init_states["DXt"] = 0.0
    init_states["DXQsto1"] = params["pDX"]
    init_states["DXQsto2"] = 0.0
    init_states["DXQgut"] = 0.0


    init_states
end # gen_init_state
# end # Dalla_Man module
