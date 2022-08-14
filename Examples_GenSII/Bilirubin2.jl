using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1+k12*x2+k13*x3+k14*x4,
  x2'(t) = k21*x1-k12*x2,
  x3'(t) = k31*x1-k13*x3,
  x4'(t) = k41*x1-k14*x4,
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

