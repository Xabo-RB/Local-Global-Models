using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = - a1*(x1(t) - x2(t)) - ka*vm*x1(t)*x5(t),
  x2'(t) =  a2*(x1(t) - x2(t)),
  x3'(t) = - b1*(x3(t) - x4(t)) - kc*vm*x3(t)*x5(t),
  x4'(t) = b2*(x3(t) - x4(t)),
  x5'(t) = ka*x5(t)^2*(a1*(x1(t) - x2(t)) + ka*vm*x1(t)*x(t)) + kc*x5(t)^2*(b1*(x3(t) - x4(t)) + kc*vm*x3(t)*x5(t)),
  y1(t) = x1(t),
  y2(t) = x4(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
