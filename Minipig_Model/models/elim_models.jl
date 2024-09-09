@parameters DXDose BW t_peak k_1 k_2 DXk_e1 DXk_e4 β_x k_x
@variables t DX_l(t) DX_p(t)

function setconstprm(symprm,dictprm,ipt)
    """
        Sets list of parameters as constant from a quote of parameters initialization
    """
    for sprm in symprm
        :($(sprm) = $(dictprm)[Symbol($(sprm))][$(ipt)]) |> eval 
    end
end

setconstprm([:BW,:t_peak,:DXDose],all_prm,1)

inj_kinetics() = IfElse.ifelse(t<t_peak, DXDose*t/(t_peak*BW),0.0)
ela_kinetics(DXp0) = DXp0*β_x*k_x^β_x*t^(β_x-1)*exp(-(k_x*t)^β_x)

# two compartments model with elimination in plasma (1st order, DXke1) and remote compartment (0 order, DXke4) 
crn_2comp_ke4 = @reaction_network begin
    inj_kinetics(), ∅ --> DX_p
    k_1, DX_p --> DX_l 
    k_2, DX_l --> DX_p
    DXk_e4, DX_l --> ∅
    DXk_e1*DX_p, DX_p ⇒ ∅
end k_1 k_2 DXk_e1 DXk_e4

# two compartments model without elimination in remote compartment (ie. without DXke4)
crn_2comp = @reaction_network begin
    inj_kinetics(), ∅ --> DX_p
    k_1, DX_p --> DX_l 
    k_2, DX_l --> DX_p
    DXk_e1*DX_p, DX_p ⇒ ∅
end k_1 k_2 DXk_e1

# two compartments model with 1st order elimination in remote compartment and no elimination in vascular compartment
crn_2comp_remote_elim = @reaction_network begin
    inj_kinetics(), ∅ --> DX_p
    k_1, DX_p --> DX_l 
    k_2, DX_l --> DX_p
    DXk_e1*DX_l, DX_l ⇒ ∅
end k_1 k_2 DXk_e1

# two compartments model with 1st order elimination in remote compartment and no elimination in vascular compartment
crn_2comp_sans_inj = @reaction_network begin
        k_1, DX_p --> DX_l 
        k_2, DX_l --> DX_p
        DXk_e1*DX_l, DX_l ⇒ ∅
end k_1 k_2 DXk_e1

# power-exponetial model
function gen_crn_powexp(DXp0)  
    @reaction_network begin
        ela_kinetics($DXp0), DX_p ⇒ ∅
    end β_x k_x
end 