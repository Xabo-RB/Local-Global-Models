using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) =  - x1(t)*(p1 + p2/(x1(t)+p3)) + p4*x2(t) + p5*u(t),
    x2'(t) = p1*x1(t) - (p4+p6)*x2(t),
    y1(t) = p7*x1(t)
)

@time println(assess_identifiability(ode))