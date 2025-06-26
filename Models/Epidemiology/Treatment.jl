using SIAN, Logging

ode = @ODEmodel(
  S'(t) = -b * S(t) * In(t) / N(t) - d * b * S(t) * Tr(t) / N(t),
  In'(t) = b * S(t) * In(t) / N(t) + d * b * S(t) * Tr(t) / N(t) - (a + g) * In(t),
  Tr'(t) = g * In(t) - nu * Tr(t),
  N'(t) = 0,
  y1(t) = Tr(t),
  y2(t) = N(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3))

#N como una entrada u(t)
using SIAN, Logging

ode = @ODEmodel(
  S'(t) = -b * S(t) * In(t) / u(t) - d * b * S(t) * Tr(t) / u(t),
  In'(t) = b * S(t) * In(t) / u(t) + d * b * S(t) * Tr(t) / u(t) - (a + g) * In(t),
  Tr'(t) = g * In(t) - nu * Tr(t),
  y1(t) = Tr(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
