using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -beta*I(t)*(S(t)/N(t)),
    E'(t) =  beta*I(t)*(S(t)/N(t)) - alpha*E(t),
    I'(t) = alpha*E(t) - lambda*I(t),
    R'(t) = lambda*I(t),
    y1(t) = I(t)
)


println(find_ioequations(ode))


using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -b * S(t) * In(t) / N(t),
    E'(t) = b * S(t) * In(t) / N(t) - nu * E(t),
    In'(t) = nu * E(t) - a * In(t),
    N'(t) = 0,
    y1(t) = In(t),
    y2(t) = N(t)
)


println(find_ioequations(ode))