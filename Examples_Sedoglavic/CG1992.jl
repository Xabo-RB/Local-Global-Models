using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = x1(t)*(p1*x2(t)/(p2+x2(t)) - p3) + p4*u(t),
    x2'(t) = - p1*x1(t)*x2(t)/(p5*(p2+x2(t))) + p6*u(t),
    y1(t) = x1(t)
)

@time println(assess_identifiability(ode))