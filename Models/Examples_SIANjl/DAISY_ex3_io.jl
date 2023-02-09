using StructuralIdentifiability

ode = @ODEmodel(
  x1'(t) = -1 * p1 * x1(t) + x2(t) + u0(t),
  x2'(t) = p3 * x1(t) - p4 * x2(t) + x3(t),
  x3'(t) = p6 * x1(t) - p7 * x3(t),
  u0'(t) = 1,
  y(t) = x1(t),
  y2(t) = u0(t)
)


println(find_ioequations(ode))
