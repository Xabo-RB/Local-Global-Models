using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t) + u0(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  u0'(t) = 1,
  y(t) = x1(t),
  y2(t) = u0(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#PRUEBAS
# u0 entrada normal->mismo resultado
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t) + u(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  y(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t), #sin entrada
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  y(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#cambiar la localización de u -> p6 y p7 SGI, p1 p3 p4 NI
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t) + u(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  y(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#cambiar la localización de u -> ninguno identificable
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t) + u(t),
  y(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
