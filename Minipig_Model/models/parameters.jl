qall_prm = :(Dict(
        # GK
        :V_G  => [1.88 , 1.49 ],
        :k_1  => [0.065, 0.042],
        :k_2  => [0.079, 0.071],
        # IK 07
        :V_I  => [0.05  , 0.04  ],
        :m_1  => [0.190 , 0.379 ],
        :m_5  => [0.0304, 0.0526],
        :m_6  => [0.6471, 0.8118],
        :HE_b => [0.6   , 0.6   ],
        # RA 
        :k_max => [0.0558, 0.0465],
        :k_min => [0.0080, 0.0076],
        :k_abs => [0.057 , 0.023 ],
        :f     => [0.90  , 0.90  ],
        :a     => [0.82  , 0.68  ],
        :b     => [0.010 , 0.09  ],
        :BW    => [120.0 , 120.0 ],
        :Dose  => [90000.0,90000.0],
        :k_t   => [0.1,0.1],
        # EGP
        :k_p1 => [2.70  , 3.09  ],
        :k_p2 => [0.0021, 0.0007],
        :k_p3 => [0.009 , 0.005 ],
        :k_p4 => [0.0618, 0.0786],
        :k_i  => [0.0079, 0.0066],
        # U 
        :F_cns => [1.0   , 1.0   ],
        :V_m0  => [2.50  , 4.65  ],
        :V_mx  => [0.047 , 0.034 ],
        :K_m0  => [225.59, 466.21],
        :p_2U  => [0.0331, 0.0840],
        # S
        :K   => [2.30 , 0.99 ],
        :α_Y => [0.050, 0.013],
        :β_Y => [0.11 , 0.05 ],
        :γ_Y => [0.5  , 0.5  ],
        :h   => [90.0,160.0],
        # RE
        :k_e1 => [0.0005, 0.0007],
        :k_e2 => [339.0 , 269.0 ],
        # DXRE
        :k_inj  => [0.1,0.1], # for any intra-veinous injection 
        :DXDose => [30000.0, 30000.0],
        :DXk_e1 => [0.0005 , 0.0007 ],
        :DXk_e2 => [0.000  , 0.000  ],
        :DXk_e3 => [0.0005   , 0.0007   ], # liver clearance
        :DXk_e4 => [0.0005   , 0.0007   ], # delayed liver clearance 
        :β_x    => [1.0,10.0],  # power exp elimination of D-xylose
        :k_x    => [0.1 , 1.0 ], # power exp elimination of D-xylose
        :t_peak => [5.0,5.0],
        # ParaC
        :ParaCDose => [1500.0, 1500.0],
        :k_pc12    => [1.0,1.0],
        :k_pc21    => [1.0,1.0],
        :k_pcel    => [0.1,0.1],
        :V_PC      => [5.0,5.0],
        # ELASHOFF
        :β_ela => [1.03 ,1.03],
        :k_ela => [0.015,1.0],
        # Basal values
        :G_b   => [90.0,160.0],
        :EGP_b => [1.92, 2.01],
        :I_b   => [25.0,50.0],
        # duplicated params for simultaneous estimations
        :k_abs1 => [0.057 , 0.023 ],
        :k_abs2 => [0.057 , 0.023 ],
    )
)

all_prm = eval(qall_prm)

# for the sake of all_u0
fS_b(m_6,HE_b,m_5) = (m_6 - HE_b)/m_5
all_prm[:S_b] = fS_b.(all_prm[:m_6],all_prm[:HE_b],all_prm[:m_5])
fm_4(S_b,I_b,V_I,HE_b) = (2.0/5.0)*(S_b/(I_b*V_I))*(1.0 - HE_b)
fm_2(S_b,I_b,V_I,HE_b) = (S_b/(I_b*V_I)-(fm_4(S_b,I_b,V_I,HE_b)/(1.0 - HE_b)))*(1.0 - HE_b)/HE_b

addu0 = Dict(:xyz => [0.0,0.0],:xyz2 => [0.0,0.0]) # for salinari or any rn that adds new variables
qall_u0 = :(Dict(
        # GK
        :G_p => all_prm[:G_b].*all_prm[:V_G],
        :G_t => (all_prm[:F_cns].-all_prm[:EGP_b].+all_prm[:k_1].*(all_prm[:G_b].*all_prm[:V_G]))./all_prm[:k_2],
        # IK
        :I_p => all_prm[:I_b].* all_prm[:V_I],
        :I_l => (all_prm[:I_b].* all_prm[:V_I]).*(fm_4.(fS_b.(all_prm[:m_6],all_prm[:HE_b],all_prm[:m_5]),all_prm[:I_b],all_prm[:V_I],all_prm[:HE_b]).+fm_2.(fS_b.(all_prm[:m_6],all_prm[:HE_b],all_prm[:m_5]),all_prm[:I_b],all_prm[:V_I],all_prm[:HE_b]))./all_prm[:m_1],
        # RA
        :Q_sto1 => all_prm[:Dose],
        :Q_sto2 => [0.0,0.0],
        :Q_gut => [0.0,0.0],
        # EGP
        :I_1 => all_prm[:I_b],
        :I_d => all_prm[:I_b],
        :X   => [0.0,0.0],
        # S
        :Y    => [0.0,0.0],
        :I_po => fS_b.(all_prm[:m_6],all_prm[:HE_b],all_prm[:m_5])./all_prm[:γ_Y],
        # DX
        :DX_sto1 => all_prm[:DXDose],
        :DX_sto2 => [0.0,0.0],
        :DX_gut  => [0.0,0.0],
        :DX_p    => [0.0,0.0],
        :DX_t    => [0.0,0.0],
        # DX simultaneous estimation
        :DX_p1 => [0.0,0.0],
        :DX_p2 => [0.0,0.0],
        :DX_t1 => [0.0,0.0],
        :DX_t2 => [0.0,0.0],
        :DX_l1 => [0.0,0.0],
        :DX_l2 => [0.0,0.0],
        :DX_gut2 => [0.0,0.0],
        :DX_gut1 => [0.0,0.0],
        :DX_sto12 => all_prm[:DXDose],
        :DX_sto11 => all_prm[:DXDose],
        :DX_sto21 => [0.0,0.0],
        :DX_sto22 => [0.0,0.0],
        :DX_inj   => all_prm[:DXDose],#mg
        # ParaC
        :ParaC_sto1 => all_prm[:ParaCDose],
        :ParaC_gut  => [0.0,0.0],
        :ParaC_p    => [0.0,0.0],
        :ParaC_t    => [0.0,0.0],
        # ParaC simultaneous estimations
        :ParaC_gut1   => [0.0,0.0],
        :ParaC_gut2   => all_prm[:ParaCDose],
        :ParaC_p1     => [0.0,0.0],
        :ParaC_t1     => [0.0,0.0],
        :ParaC_p2     => [0.0,0.0],
        :ParaC_t2     => [0.0,0.0],
        $addu0...
    )
)

all_u0 = eval(qall_u0)
   
function gen_u0(all_prm,symprm,qall_u0)
    """
        Generates the initial state with:
         - the dictionary of parameters 
         - the parameters names as Symbols
    """
    newu0 = eval(qall_u0)
    newu0
end




function update_u0(all_prm,symprm,symu0)
    """
        Generates the initial state with:
         - the dictionary of parameters 
         - the parameters names as Symbols
    """ 
    newu0 =  Dict(
        # GK
        :G_p => all_prm[:G_b].*all_prm[:V_G],
        :G_t => (all_prm[:F_cns].-all_prm[:EGP_b].+all_prm[:k_1].*(all_prm[:G_b].*all_prm[:V_G]))./all_prm[:k_2],
        # IK
        :I_p => all_prm[:I_b].* all_prm[:V_I],
        # RA
        :Q_sto1 => all_prm[:Dose],
        :Q_sto2 => [0.0,0.0],
        :Q_gut => [0.0,0.0],
        # EGP
        :I_1 => all_prm[:I_b],
        :I_d => all_prm[:I_b],
        :X   => [0.0,0.0],
        # S
        :Y    => [0.0,0.0],
        # DX
        :DX_sto1 => all_prm[:DXDose],
        :DX_sto2 => [0.0,0.0],
        :DX_gut  => [0.0,0.0],
        :DX_p    => [0.0,0.0],
        :DX_t    => [0.0,0.0]
    )
    if(:I_po in symu0)
        newu0[:I_po] = fS_b.(all_prm[:m_6],all_prm[:HE_b],all_prm[:m_5])./all_prm[:γ_Y]
        newu0[:I_l] = (all_prm[:I_b].* all_prm[:V_I]).*(fm_4.(fS_b.(all_prm[:m_6],all_prm[:HE_b],all_prm[:m_5]),all_prm[:I_b],all_prm[:V_I],all_prm[:HE_b]).+fm_2.(fS_b.(all_prm[:m_6],all_prm[:HE_b],all_prm[:m_5]),all_prm[:I_b],all_prm[:V_I],all_prm[:HE_b]))./all_prm[:m_1]
        #else 
    end
newu0
end



;