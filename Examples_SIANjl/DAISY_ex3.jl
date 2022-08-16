#MODELO IDENTIFICABLE, SIN P2
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t) + u0(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  u0'(t) = 1,
  y(t) = x1(t),
  y2(t) = u0(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#MODELO NO IDENTIFICABLE, CON P2 y P5
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + p2*x2(t) + u0(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + p5*x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  u0'(t) = 1,
  y(t) = x1(t),
  y2(t) = u0(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)


#MODELO CON P2 Y DEFINIENDO U COMO UNA ENTRADA
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + p2*X2(t) + u(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + p5*x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  y(t) = x1(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

#MODELO CON MENOS PARAMETROS CON U0 = CONSTANTE -> MISMO RESULTADOS
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t) + u0(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  u0'(t) = 1,
  y(t) = x1(t),
  y2(t) = u0(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

#MODELO CON MENOS PARAMETROS CON U0 = CONSTANTE -> MISMO RESULTADOS
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t) + u0(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  u0'(t) = 0,
  y(t) = x1(t),
  y2(t) = u0(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

#MODELO SIN ENTRADA
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + p2*x2(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + p5*x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  y(t) = x1(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

