using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = k3*x3(t) - r3*x1(t) - k1*x1(t)*(T0-x2(t)) + r1*x2(t),
    x2'(t) = k1*x1(t)*(T0-x2(t)) - (r1+k2)*x2(t),
    x3'(t) = r3*x1(t) - (k3+k4)*x3(t) + k2*x2(t),
    y1(t) = k*(x2(t) + x3(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
