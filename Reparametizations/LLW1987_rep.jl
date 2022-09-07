using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + u(t),
    x2'(t) = - p3*x2(t) + p4*u(t),
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+x2(t))*u(t),
    y1(t) = x3(t)
)

@time println(assess_identifiability(ode))
