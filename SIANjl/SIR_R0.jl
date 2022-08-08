using SIAN, Logging

ode = @ODEmodel(
  S'(t) = -b * In(t) * S(t),
  In'(t) = b * In(t) * S(t) - g * In(t),
  R'(t) = g * In(t),
  aux'(t) = 0,
  y1(t) = In(t),
  y2(t) = b // g + aux(t)
)

res = identifiability_ode(ode, [aux]; p = 0.99, p_mod = 0, nthrds = 1)

println(res)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))