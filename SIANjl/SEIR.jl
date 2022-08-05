using SIAN, Logging
@info "Setting up the problem"

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