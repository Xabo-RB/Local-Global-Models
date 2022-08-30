using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = t1*x1(t)^2+t2*x1(t)*x2(t)+u(t),
    x2'(t) = t3*x1(t)^2+t4*x1(t)*x2(t),
    y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
