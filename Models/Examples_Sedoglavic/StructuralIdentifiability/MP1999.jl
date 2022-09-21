using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = x3(t) - x2(t)*u(t),
    x2'(t) = u(t) - x2(t),
    x3'(t) = x2(t) - x1(t) + 2*x2(t)*(u(t)-x2(t)),
    y1(t) = x1(t)+(x2(t)^2)/2
)

@time println(assess_identifiability(ode))
