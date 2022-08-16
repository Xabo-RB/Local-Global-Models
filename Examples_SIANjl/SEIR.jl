using SIAN, Logging

ode = @ODEmodel(
  S'(t) = -b * S(t) * In(t) / N(t),
  E'(t) = b * S(t) * In(t) / N(t) - nu * E(t),
  In'(t) = nu * E(t) - a * In(t),
  N'(t) = 0,
  y1(t) = In(t),
  y2(t) = N(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#PRUEBA CON INPUT -> mismo resultado
using SIAN, Logging

ode = @ODEmodel(
  S'(t) = -b * S(t) * In(t) / u(t),
  E'(t) = b * S(t) * In(t) / u(t) - nu * E(t),
  In'(t) = nu * E(t) - a * In(t),
  y1(t) = In(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

#DEFINIR LA CONSTANTE CONOCIDA N como un parámetro y no como un estado. -> mismo resultado
using SIAN, Logging

ode = @ODEmodel(
  S'(t) = -b * S(t) * In(t) / N,
  E'(t) = b * S(t) * In(t) / N - nu * E(t),
  In'(t) = nu * E(t) - a * In(t),
  y1(t) = In(t),
  y2(t) = N
)

  res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)
  
  println(res)


#SIN N(T) DE NINGUNA FORMA  
using SIAN, Logging

ode = @ODEmodel(
  S'(t) = -b * S(t) * In(t),
  E'(t) = b * S(t) * In(t) - nu * E(t),
  In'(t) = nu * E(t) - a * In(t),
  y1(t) = In(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)
