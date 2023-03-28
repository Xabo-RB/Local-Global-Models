using StructuralIdentifiability

ode = @ODEmodel(
  S'(t) = -b * S(t) * In(t) / N(t) - d * b * S(t) * Tr(t) / N(t),
  In'(t) = b * S(t) * In(t) / N(t) + d * b * S(t) * Tr(t) / N(t) - (a + g) * In(t),
  Tr'(t) = g * In(t) - nu * Tr(t),
  N'(t) = 0,
  y1(t) = Tr(t),
  y2(t) = N(t)
)

println(find_ioequations(ode))

