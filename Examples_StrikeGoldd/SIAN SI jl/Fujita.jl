using SIAN, Logging

ode = @ODEmodel(
    EGFR'(t) = 68190.0*p1 - EGF_u(t)*EGFR(t)*p2 + EGF_EGFR(t)*p3 - EGFR(t)*p1,
    pEGFR'(t) = EGF_EGFR(t)*p4 - pEGFR(t)*p5 + pEGFR_Akt(t)*p6 + pEGFR_Akt(t)*p7 - Akt(t)*pEGFR(t)*p8,
    pEGFR_Akt'(t) = Akt(t)*pEGFR(t)*p8 - pEGFR_Akt(t)*p7 - pEGFR_Akt(t)*p6,
    Akt'(t) = pAkt(t)*p9 + pEGFR_Akt(t)*p6 - Akt(t)*pEGFR(t)*p8,
    pAkt'(t) = pAkt_S6(t)*p10 - pAkt(t)*p9 + pAkt_S6(t)*p12 + pEGFR_Akt(t)*p7 - S6(t)*pAkt(t)*p13,
    S6'(t) = pAkt_S6(t)*p10 + pS6(t)*p11 - S6(t)*pAkt(t)*p13,
    pAkt_S6'(t) = S6(t)*pAkt(t)*p13 - pAkt_S6(t)*p12 - pAkt_S6(t)*p10,
    pS6'(t) = pAkt_S6(t)*p12 - pS6(t)*p11,
    EGF_EGFR'(t) = EGF_u(t)*EGFR(t)*p2 - EGF_EGFR(t)*p3 - EGF_EGFR(t)*p4,
    y1(t) = p14*(pEGFR(t) + pEGFR_Akt(t)),
    y2(t) = p15*(pAkt(t) + pAkt_S6(t)),
    y3(t) = pS6(t)*p16
)

#p1 = EGFR_turnover
#p2 = reaction_1_k1
#p3 = reaction_1_k2
#p4 = reaction_9_k1
#p5 = reaction_4_k1
#p6 = reaction_2_k2

#p9 = reaction_7_k1
#p10 = reaction_5_k2
#p11 = reaction_8_k1
#p12 = reaction_6_k1
#p13 = reaction_5_k1
#p14 = scaleFactor_pEGFR
#p15 = scaleFactor_pAkt
#p16 = scaleFactor_pS6




@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
