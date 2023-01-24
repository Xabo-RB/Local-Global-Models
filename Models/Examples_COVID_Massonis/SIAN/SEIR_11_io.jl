# https://github.com/ryansmcgee/seirsplus#model



using StructuralIdentifiability


ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*S(t)*I(t)/N(t)+xi*R(t),
    E'(t) = beta*S(t)*I(t)/N(t)-sigma*E(t),
    I'(t) = sigma*E(t)-gamma*I(t)-mu*I(t),
    R'(t) = gamma*I(t)-xi*R(t),
    F'(t) = mu*I(t),
    y1(t) = F(t),
    y2(t) = N(t)
)

#@time println(assess_identifiability(ode))

println(find_ioequations(ode))