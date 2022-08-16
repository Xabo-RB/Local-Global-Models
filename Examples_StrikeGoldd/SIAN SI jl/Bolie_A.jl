using SIAN, Logging
ode = @ODEmodel(
    q1'(t) = p1*q1(t)-p2*q2(t) + u(t),
    q2'(t) = p3*q2(t)+p4*q1(t),
    y1(t) = q1(t)/Vp
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
