using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = u2(t)+u1(t)-p2*x1(t)-p1*x2(t),
    x2'(t) = x1(t)-x10,
    y1(t) = x1(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))