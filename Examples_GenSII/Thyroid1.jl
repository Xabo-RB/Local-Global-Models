using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = k13*x3(t)+k12*x2(t)-(k21+k31)*x1(t) + u(t),
  x2'(t) =  k21*x1(t)-(k12+k02)*x2(t),
  x3'(t) = k31*x1(t)-(k13+k03)*x3(t),
  y1(t) = x1(t)/V1
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
