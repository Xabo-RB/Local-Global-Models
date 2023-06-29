using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = -b*x1(t) + a/(AA+k*x2(t)),
    x2'(t) = gamma*x3(t) - delta*x2(t),
    x3'(t) = u(t),
    y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
