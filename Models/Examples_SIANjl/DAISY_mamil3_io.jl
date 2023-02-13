using StructuralIdentifiability

ode = @ODEmodel(
  x1'(t) = -(a21 + a31 + a01) * x1(t) + a12 * x2(t) + a13 * x3(t) + u(t),
  x2'(t) = a21 * x1(t) - a12 * x2(t),
  x3'(t) = a31 * x1(t) - a13 * x3(t),
  y(t) = x1(t)
)


println(find_ioequations(ode))
