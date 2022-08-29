using SIAN, Logging

ode = @ODEmodel(
    EGFR'(t) = 68190.0*EGFR_turnover - EGF_u(t)*EGFR(t)*reaction_1_k1 + EGF_EGFR(t)*reaction_1_k2 - EGFR(t)*EGFR_turnover,
    pEGFR'(t) = EGF_EGFR(t)*reaction_9_k1 - pEGFR(t)*reaction_4_k1 + pEGFR_Akt(t)*reaction_2_k2 + pEGFR_Akt(t)*reaction_3_k1 - Akt(t)*pEGFR(t)*reaction_2_k1,
    pEGFR_Akt'(t) = Akt(t)*pEGFR(t)*reaction_2_k1 - pEGFR_Akt(t)*reaction_3_k1 - pEGFR_Akt(t)*reaction_2_k2,
    Akt'(t) = pAkt(t)*reaction_7_k1 + pEGFR_Akt(t)*reaction_2_k2 - Akt(t)*pEGFR(t)*reaction_2_k1,
    pAkt'(t) = pAkt_S6(t)*reaction_5_k2 - pAkt(t)*reaction_7_k1 + pAkt_S6(t)*reaction_6_k1 + pEGFR_Akt(t)*reaction_3_k1 - S6(t)*pAkt(t)*reaction_5_k1,
    S6'(t) = pAkt_S6(t)*reaction_5_k2 + pS6(t)*reaction_8_k1 - S6(t)*pAkt(t)*reaction_5_k1,
    pAkt_S6'(t) = S6(t)*pAkt(t)*reaction_5_k1 - pAkt_S6(t)*reaction_6_k1 - pAkt_S6(t)*reaction_5_k2,
    pS6'(t) = pAkt_S6(t)*reaction_6_k1 - pS6(t)*reaction_8_k1,
    EGF_EGFR'(t) = EGF_u(t)*EGFR(t)*reaction_1_k1 - EGF_EGFR(t)*reaction_1_k2 - EGF_EGFR(t)*reaction_9_k1,
    y1(t) = scaleFactor_pEGFR*(pEGFR(t) + pEGFR_Akt(t)),
    y2(t) = scaleFactor_pAkt*(pAkt(t) + pAkt_S6(t));
    y3(t) = pS6(t)*scaleFactor_pS6
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
