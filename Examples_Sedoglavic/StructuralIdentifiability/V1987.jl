using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -x1(t)*(k1+k2*x4(t)) + k5*x3(t)*x4(t),
    x2'(t) = k2*x1(t)*x4(t) - (k3+k4)*x2(t),
    x3'(t) = k4*x2(t) - k5*x3(t)*x4(t),
    x4'(t) = x1(t)*(k1+k2*x4(t)) +2*k3*x2(t) - k5*x3(t)*x4(t),
    y1(t) = x1(t),
    y2(t) = x2(t)
)

@time println(assess_identifiability(ode))
