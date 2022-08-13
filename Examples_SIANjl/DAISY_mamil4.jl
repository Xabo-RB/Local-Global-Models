
using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = -k01 * x1(t) + k12 * x2(t) + k13 * x3(t) + k14 * x4(t) - k21 * x1(t) - k31 * x1(t) - k41 * x1(t) + u(t),
    x2'(t) = -k12 * x2(t) + k21 * x1(t),
    x3'(t) = -k13 * x3(t) + k31 * x1(t),
    x4'(t) = -k14 * x4(t) + k41 * x1(t),
    y(t) = x1(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, nthrds = 1)

println(res)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#SIN ENTRADA
using SIAN, Logging

ode = @ODEmodel(
    x1'(t) = -k01 * x1(t) + k12 * x2(t) + k13 * x3(t) + k14 * x4(t) - k21 * x1(t) - k31 * x1(t) - k41 * x1(t),
    x2'(t) = -k12 * x2(t) + k21 * x1(t),
    x3'(t) = -k13 * x3(t) + k31 * x1(t),
    x4'(t) = -k14 * x4(t) + k41 * x1(t),
    y(t) = x1(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, nthrds = 1)

println(res)