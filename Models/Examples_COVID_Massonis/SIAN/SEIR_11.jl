# https://github.com/ryansmcgee/seirsplus#model

using SIAN, Logging

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

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#PRUEBAS
#sigue igual
using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t)+xi*R(t),
    E'(t) = beta*S(t)*I(t)-sigma*E(t),
    I'(t) = sigma*E(t)-gamma*I(t)-mu*I(t),
    R'(t) = gamma*I(t)-xi*R(t),
    F'(t) = mu*I(t),
    y1(t) = F(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#TODOS LOS SLI A NI
using SIAN, Logging

ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*S(t)*I(t)/N(t)+xi*R(t),
    E'(t) = beta*S(t)*I(t)/N(t)-E(t),#ELIMINÉ SIGMA
    I'(t) = sigma*E(t)-gamma*I(t)-mu*I(t),
    R'(t) = gamma*I(t)-xi*R(t),
    F'(t) = mu*I(t),
    y1(t) = F(t),
    y2(t) = N(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#La mitad a NI y la mitad a sgi
using SIAN, Logging

ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*S(t)*I(t)/N(t)+xi*R(t),
    E'(t) = beta*S(t)*I(t)/N(t)-sigma*E(t),
    I'(t) = sigma*E(t)-gamma*I(t)-I(t), #elimino MU
    R'(t) = gamma*I(t)-xi*R(t),
    F'(t) = mu*I(t),
    y1(t) = F(t),
    y2(t) = N(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
