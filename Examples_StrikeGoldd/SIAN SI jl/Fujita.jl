using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = 68190.0*p1 - u(t)*x1(t)*p2 + x9(t)*p3 - x1(t)*p1,
    x2'(t) = x9(t)*p4 - x2(t)*p5 + x3(t)*p6 + x3(t)*p7 - x4(t)*x2(t)*p8,
    x3'(t) = x4(t)*x2(t)*p8 - x3(t)*p7 - x3(t)*p6,
    x4'(t) = x5(t)*p9 + x3(t)*p6 - x4(t)*x2(t)*p8,
    x5'(t) = x6(t)*p10 - x5(t)*p9 + x6(t)*p12 + x3(t)*p7 - x7(t)*x5(t)*p13,
    x7'(t) = x6(t)*p10 + x8(t)*p11 - x7(t)*x5(t)*p13,
    x6'(t) = x7(t)*x5(t)*p13 - x6(t)*p12 - x6(t)*p10,
    x8'(t) = x6(t)*p12 - x8(t)*p11,
    x9'(t) = u(t)*x1(t)*p2 - x9(t)*p3 - x9(t)*p4,
    y1(t) = p14*(x2(t) + x3(t)),
    y2(t) = p15*(x5(t) + x6(t)),
    y3(t) = x8(t)*p16
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

#x1 = EGFR
#x2 = pEGFR
#x3 = pEGFR_Akt
#x4 = Akt
#x5 = pAkt
#x6 = pAkt_S6
#x7 = S6
#x8 = pS6
#x9 = EGF_EGFR


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
