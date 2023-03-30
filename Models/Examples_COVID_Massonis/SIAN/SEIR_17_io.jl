#Healthcare impact of COVID-19 epidemic in India: A stochastic mathematical model

using StructuralIdentifiability
ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t),
    E'(t) = beta*S(t)*I(t)-epsilon*E(t),
    I'(t) = epsilon*E(t)-gamma*I(t)-d*I(t)-q*I(t),
    Q'(t) = q*I(t)-qt*Q(t)-d*Q(t),
    R'(t) = gamma*I(t)+qt*Q(t),
    D'(t) = d*I(t)+d*Q(t),
    y1(t) = D(t)
)
println(find_ioequations(ode))
