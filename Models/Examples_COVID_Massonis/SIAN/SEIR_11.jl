# https://github.com/ryansmcgee/seirsplus#model

using SIAN, Logging

ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*S(t)*I(t)/N(t)+xi*R(t),
    E'(t) = beta*S(t)*I(t)/N(t)-sigma*E(t),
    I'(t) = sigma*E(t)-gamma*I(t)-mu*I(t),
    R'(t) = gamma*I(t)-xi*R(t),
    F'(t) = mu*I(t),
    y1(t) = F(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
