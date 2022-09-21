using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = u(t) - (c1+c2)*x1(t),
    x2'(t) = c1*x1(t)-(c3+c6+c7)*x2(t) + c5*x4(t),
    x3'(t) = c2*x1(t) + c3*x2(t) - c4*x3(t),
    x4'(t) = c6*x2(t) - c5*x4(t),
    y1(t) = c8*x3(t),
    y2(t) = c9*x2(t)
)

@time println(assess_identifiability(ode))