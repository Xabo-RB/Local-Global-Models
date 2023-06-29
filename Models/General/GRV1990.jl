using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = p1*x1(t)^2 + p2*x1(t)*x2(t) + u(t),
    x2'(t) = p3*x1(t)^2 + p4*x1(t)*x2(t),
    y1(t) = x1(t)
)

@time println(assess_identifiability(ode))
