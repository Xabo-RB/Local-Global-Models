using StructuralIdentifiability

ode = @ODEmodel(
    x4'(t) = -k5 * x4(t) // (k6 + x4(t)),
    x5'(t) = k5 * x4(t) // (k6 + x4(t)) - k7 * x5(t) / (k8 + x5(t) + x6(t)),
    x6'(t) = k7 * x5(t) // (k8 + x5(t) + x6(t)) - k9 * x6(t) * (k10 - x6(t)) // k10,
    x7'(t) = k9 * x6(t) * (k10 - x6(t)) // k10,
    y1(t) = x4(t),
    y2(t) = x5(t)
)

println(find_ioequations(ode))
