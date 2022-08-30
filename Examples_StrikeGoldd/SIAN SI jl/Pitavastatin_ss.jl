using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = k3*x3(t) - r3*x1(t) - Vm*x1(t)/(Km+x1(t)),
    x3'(t) = r3*x1(t) - (k3+k4)*x3(t) + Vm*x1(t)/(Km+x1(t)),
    y1(t) = k*( T0*x1(t)/(Km+x1(t)) + x3(t))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
