using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k1*x1(t))/(kM+x1(t))*u1(t),
  x2'(t) = (k1*x1(t))/(kM+x1(t))*u1(t) - (k2*x2(t))/(kM+x2(t))*u2(t),
  x3'(t) = (k2*x2(t))/(kM+x2(t))*u2(t) - (k3*x3(t))/(kM+x3(t))*u3(t),
  x4'(t) = (k2*x2(t))/(kM+x2(t))*u2(t) + (k3*x3(t))/(kM+x3(t))*u3(t) - (k4*x4(t))/(kM+x4(t))*u4(t),
  x5'(t) = (k4*x4(t))/(kM+x4(t))*u4(t),
  y1(t) = x1(t),
  y2(t) = x2(t),
  y3(t) = x3(t),
  y4(t) = x4(t),
  y5(t) = x5(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
