using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = x1(t)*(a-b*x1(t))-c*x1(t),
    y1(t) = d*x1(t)
)

@time println(assess_identifiability(ode))
