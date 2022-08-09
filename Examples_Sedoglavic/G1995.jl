using StructuralIdentifiability

ode = @ODEmodel(
    M'(t) = vs*KI^4/(KI^4 + PN(t)^4) - (vm*M(t)/(Km + M(t))),
    P0'(t) = (ks*M(t) - V1*P0(t)/(K1+P0(t)) + V2*P1(t)/(K2+P1(t))) ,
    P1'(t) = (V1*P0(t)/(K1+P0(t)) - P1(t)*(V2/(K2+P1(t)) + V3/(K3+P1(t))) + V4*P2(t)/(K4+P2(t))),
    P2'(t) = (V3*P1(t)/(K3+P1(t)) - P2(t)*(V4/(K4+P2(t)) + k1 + vd/(Kd+P2(t))) + k2*PN(t)),
    PN'(t) = k1*P2(t) - k2*PN(t),
    y1(t) = PN(t)
)

@time println(assess_identifiability(ode))