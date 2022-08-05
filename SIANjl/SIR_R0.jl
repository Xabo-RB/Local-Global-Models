using SIAN, Logging
@info "Setting up the problem"

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