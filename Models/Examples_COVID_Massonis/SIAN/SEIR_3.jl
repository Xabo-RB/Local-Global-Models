# Liangrong Peng et al. (2020) Epidemic analysis of COVID-19 in 
# China by dynamical modeling 


using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) =  beta*I(t)*S(t)/N(t)-gamma*E(t),
    I'(t) = gamma*E(t)-delta*I(t),
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#PRUEBAS PARA ANÁLISIS
#Mismos resultados sin N
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = beta*I(t)*S(t)-alpha*S(t),
    E'(t) =  beta*I(t)*S(t)-gamma*E(t),
    I'(t) = gamma*E(t)-delta*I(t),
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#No funciona
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = beta*I(t)*S(t)/N-alpha*S(t),
    E'(t) =  beta*I(t)*S(t)/N-gamma*E(t),
    I'(t) = gamma*E(t)-delta*I(t),
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#SLI todos, no son locales por el signo, lo son por la simetría
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) =  beta*I(t)*S(t)/N(t) + gamma*E(t), #signo cambiado
    I'(t) = gamma*E(t) + delta*I(t), #signo cambiado
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#TODAS GLOBALES SALVO p
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) = beta*I(t)*S(t)/N(t)-gamma + E(t), #deshice la simetría al separar las variables
    I'(t) = gamma + E(t)-delta*I(t),
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#TODOS GLOBALES SALVO P
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) = beta*I(t)*S(t)/N(t)-gamma*E(t), 
    I'(t) = gamma*E(t)-delta+I(t),#deshice la simetría al separar las variables
    R'(t) = delta+I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#TODOS A SGI O NI
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) = beta*I(t)*S(t)/N(t)-gamma*E(t), 
    I'(t) = gamma*E(t),#-delta*I(t), ELIMINÉ ESTO
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) = beta*I(t)*S(t)/N(t),#-gamma*E(t), 
    I'(t) = gamma*E(t)-delta*I(t),
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) =  beta*I(t)*S(t)/N(t)-gamma*E(t),
    I'(t) = gamma*E(t)-delta*I(t),
    R'(t) = delta*I(t)-alpha*Q(t)-k*Q(t),
    Q'(t) = alpha*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
