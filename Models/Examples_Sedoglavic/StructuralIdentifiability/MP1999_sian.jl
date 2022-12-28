using SIAN, Logging

#TODOS GLOBALES PERO NO TIENE PARÁMETROS
ode = @ODEmodel(
    x1'(t) = x3(t) - x2(t)*u(t),
    x2'(t) = u(t) - x2(t),
    x3'(t) = x2(t) - x1(t) + 2*x2(t)*(u(t)-x2(t)),
    y1(t) = x1(t)+(x2(t)*x2(t))/2
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = x3(t) - x2(t)*u(t),
    x2'(t) = u(t) - x2(t),
    x3'(t) = x2(t) - x1(t) + 2*x2(t)*(u(t)-x2(t)),
    y1(t) = x1(t)+(x2(t)^2)/2
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))