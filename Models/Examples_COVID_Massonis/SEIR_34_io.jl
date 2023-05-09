# Chintis, N. 2017. Introduction to SEIR Models


using StructuralIdentifiability


ode = @ODEmodel(
    A'(t) = 0,
    N'(t) = 0,
    S'(t) = A(t)-r*beta*S(t)*I(t)/N(t)-mu*S(t),
    E'(t) = r*beta*S(t)*I(t)/N(t)-epsilon*E(t)-mu*E(t),
    I'(t) = epsilon*E(t)-gamma*I(t)-mu*I(t),
    R'(t) = gamma*I(t)-mu*R(t),
    y1(t) = K*I(t),
    y2(t) = A(t),
    y3(t) = N(t)
)

#@time println(assess_identifiability(ode))

println(find_ioequations(ode))