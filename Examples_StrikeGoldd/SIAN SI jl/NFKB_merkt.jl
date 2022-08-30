using SIAN, Logging

ode = @ODEmodel(
    x1'(t) =  k11*x10(t)-(k1*u1(t)/(1+k0*u1(t))+k1p)*x1(t),
    x2'(t) = (k1*u1(t)/(1+k0*u1(t))+k1p)*x1(t)-k2*x2(t),
    x3'(t) =  k2*x2(t)-k3*x3(t),
    x4'(t) = k2*x2(t) - k4*x4(t),
    x5'(t) = k3*rhovol*x3(t) - k5*x5(t),
    x6'(t) = k5*x5(t) - k10*x9(t)*x6(t),
    x7'(t) = k6*x6(t) - k7*x7(t),
    x8'(t) = k8*x7(t) - k9*x8(t),
    x9'(t) = k9*rhovol*x8(t) - k10*x9(t)*x6(t),
    x10'(t) = k10*x9(t)*x6(t) - k11*rhovol*x10(t),
    y1(t) = s1*(x1(t)+x2(t)+x3(t))+I0cyt,
    y2(t) = s2*(x10(t)+x5(t)+x6(t))+I0nuc,
    y3(t) = s3*(x2(t)+x3(t)),
    y4(t) = s4*(x2(t)+x4(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
