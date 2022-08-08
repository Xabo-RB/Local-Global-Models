using SIAN, Logging

ode = @ODEmodel(
  G'(t) = -(p1 + X(t)) * G(t) + p1 * Gb(t) + v * R(t),
  X'(t) = -p2 * X(t) + p3 * (u(t) - Ib(t)),
  R'(t) = k,
  Ib'(t) = 0,
  Gb'(t) = 0,
  y1(t) = G(t),
  y2(t) = Ib(t),
  y3(t) = Gb(t)
)

res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
