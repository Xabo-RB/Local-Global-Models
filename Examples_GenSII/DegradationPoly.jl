using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = k_syn - k_deg_max*x1(t)*x2(t),
  x2'(t) = -x2(t)^2*(k_syn - k_deg_max*x1(t)*x2(t)),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))