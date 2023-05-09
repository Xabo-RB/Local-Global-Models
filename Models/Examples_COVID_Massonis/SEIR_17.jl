#Healthcare impact of COVID-19 epidemic in India: A stochastic mathematical model

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t),
    E'(t) = beta*S(t)*I(t)-epsilon*E(t),
    I'(t) = epsilon*E(t)-gamma*I(t)-d*I(t)-q*I(t),
    Q'(t) = q*I(t)-qt*Q(t)-d*Q(t),
    R'(t) = gamma*I(t)+qt*Q(t),
    D'(t) = d*I(t)+d*Q(t),
    y1(t) = D(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#PRUEBAS
#No cambia nada
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t),
    E'(t) = beta*S(t)*I(t)-epsilon*E(t),
    I'(t) = E(t)-gamma*I(t)-d*I(t)-q*I(t), #elimino epsilon
    Q'(t) = q*I(t)-qt*Q(t)-d*Q(t),
    R'(t) = gamma*I(t)+qt*Q(t),
    D'(t) = d*I(t)+d*Q(t),
    y1(t) = D(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#Se pierde el SLI
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t),
    E'(t) = beta*S(t)*I(t)-E(t),
    I'(t) = epsilon*E(t)-gamma*I(t)-d*I(t)-q*I(t),
    Q'(t) = q*I(t)-qt*Q(t)-d*Q(t),
    R'(t) = gamma*I(t)+qt*Q(t),
    D'(t) = d*I(t)+d*Q(t),
    y1(t) = D(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#No cambia nada
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t),
    E'(t) = beta*S(t)*I(t)-epsilon*E(t),
    I'(t) = epsilon*E(t)-gamma*I(t)-d*I(t)-q*I(t),
    Q'(t) = q*I(t)-qt*Q(t)-d*Q(t),
    R'(t) = I(t)+qt*Q(t), #ELIMINO GAMMA
    D'(t) = d*I(t)+d*Q(t),
    y1(t) = D(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
