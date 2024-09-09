using Reduce
eq1 = R"pk_1 * Gpb = pEGPb - pF_cns - pk_e1*Gpb - pk_e2*pk_e1 + pk_e2*Gtb"
eq2 = R"Gtb * (pV_m0- pEGPb + pF_cns - pk_e1*pk_e2 + pk_e1*Gpb) = pxK_m0*(pEGPb - pF_cns + pk_e1*pk_e2 - pk_e1*Gpb)"
sol = Algebra.solve((eq1,eq2),(:Gpb,:Gtb))

eq1 = :(pk_1 * Gpb = pEGPb - pF_cns - pk_e1*Gpb - pk_e2*pk_e1 + pk_e2*Gtb)
eq2 = :(Gtb * (pV_m0- pEGPb + pF_cns - pk_e1*pk_e2 + pk_e1*Gpb) = pxK_m0*(pEGPb - pF_cns + pk_e1*pk_e2 - pk_e1*Gpb))
sol = Algebra.solve((eq1,eq2),(:Gpb,:Gtb))

Gp = ((((pk_e1 * pk_e2 - pv_m0) * pk_1 - (pk_e2 * pxk_m0 + pv_m0) * pk_e1) - (pk_1 + 2pk_e1) * pf_cns) + (pk_1 + 2pk_e1) * pegpb + sqrt(((((((((((((((((pegpb ^ 2 * pk_1 ^ 2 - 2 * pegpb * pf_cns * pk_1 ^ 2) + 2 * pegpb * pk_1 ^ 2 * pk_e1 * pk_e2) - 2 * pegpb * pk_1 ^ 2 * pv_m0) + 4 * pegpb * pk_1 * pk_e1 ^ 2 * pk_e2 + 2 * pegpb * pk_1 * pk_e1 * pk_e2 * pxk_m0) - 2 * pegpb * pk_1 * pk_e1 * pv_m0) + pf_cns ^ 2 * pk_1 ^ 2) - 2 * pf_cns * pk_1 ^ 2 * pk_e1 * pk_e2) + 2 * pf_cns * pk_1 ^ 2 * pv_m0) - 4 * pf_cns * pk_1 * pk_e1 ^ 2 * pk_e2) - 2 * pf_cns * pk_1 * pk_e1 * pk_e2 * pxk_m0) + 2 * pf_cns * pk_1 * pk_e1 * pv_m0 + pk_1 ^ 2 * pk_e1 ^ 2 * pk_e2 ^ 2) - 2 * pk_1 ^ 2 * pk_e1 * pk_e2 * pv_m0) + pk_1 ^ 2 * pv_m0 ^ 2 + 4 * pk_1 * pk_e1 ^ 3 * pk_e2 ^ 2 + 2 * pk_1 * pk_e1 ^ 2 * pk_e2 ^ 2 * pxk_m0) - 6 * pk_1 * pk_e1 ^ 2 * pk_e2 * pv_m0) + 2 * pk_1 * pk_e1 * pk_e2 * pv_m0 * pxk_m0 + 2 * pk_1 * pk_e1 * pv_m0 ^ 2 + 4 * pk_e1 ^ 4 * pk_e2 ^ 2 + 4 * pk_e1 ^ 3 * pk_e2 ^ 2 * pxk_m0) - 4 * pk_e1 ^ 3 * pk_e2 * pv_m0) + pk_e1 ^ 2 * pk_e2 ^ 2 * pxk_m0 ^ 2 + 2 * pk_e1 ^ 2 * pk_e2 * pv_m0 * pxk_m0 + pk_e1 ^ 2 * pv_m0 ^ 2)) / (2 * (pk_1 + pk_e1) * pk_e1)

eq1 = :(pEGPb = Uidb + pF_cns)
eq2 = :(Gtb * pV_m0 = (pEGPb-pF_cns)*(pxK_m0 + Gtb))
eq3 = :(pk_1 * Gpb = pEGPb - pF_cns + pk_2*Gtb)
sol = Algebra.solve((eq1,eq2,eq3),(:Gpb,:Uidb,:Gtb))