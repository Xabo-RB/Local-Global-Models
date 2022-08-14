using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = p1*x4(t) - p4*x1(t),
    x2'(t) =  p5*x1(t) - p6*x2(t),
    x3'(t) =  p7*x2(t) - p8*x3(t),
    x4'(t) = -p3*x3(t)^(p3 - 1)*x4(t)^2*(p7*x2(t) - p8*x3(t)),
    y1(t) = x1(t),
    y2(t) = x2(t),
    y3(t) = x3(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))