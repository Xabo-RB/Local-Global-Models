using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = -x2(t) + x1(t) - x1(t)^3/3,
    x2'(t) = epsilon*alpha*x1(t) -epsilon*(x2(t)+eta) +epsilon*alpha*x1(t) ,
    y1(t) = x1(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
