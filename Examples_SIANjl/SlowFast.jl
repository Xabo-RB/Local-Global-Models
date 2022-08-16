using SIAN, Logging

ode = @ODEmodel(
  xA'(t) = -k1 * xA(t),
  xB'(t) = k1 * xA(t) - k2 * xB(t),
  xC'(t) = k2 * xB(t),
  eA'(t) = 0,
  eC'(t) = 0,
  y1(t) = xC(t),
  y2(t) = eA(t) * xA(t) + eB * xB(t) + eC(t) * xC(t),
  y3(t) = eA(t),
  y4(t) = eC(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#SIN U O CONSTANTES CONOCIDAS

using SIAN, Logging

ode = @ODEmodel(
  xA'(t) = -k1 * xA(t),
  xB'(t) = k1 * xA(t) - k2 * xB(t),
  xC'(t) = k2 * xB(t),
  y1(t) = xC(t),
  y2(t) = xA(t) + xB(t) + xC(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
