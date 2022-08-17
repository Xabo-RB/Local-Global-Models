using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = -2*x1(t)*x1(t)*k,
    y1(t) = s1*x1(t)/(1+s2*x1(t))
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
