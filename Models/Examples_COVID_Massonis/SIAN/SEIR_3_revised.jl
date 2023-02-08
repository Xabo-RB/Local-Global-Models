# Liangrong Peng et al. (2020) Epidemic analysis of COVID-19 in 
# China by dynamical modeling 

#MISMOS RESULTADOS
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) =  beta*I(t)*S(t)/N(t)-gamma*E(t),
    I'(t) = gamma*E(t)-delta*I(t),
    R'(t) = delta*I(t)-lambda*Q(t)-k*Q(t),
    Q'(t) = lambda*Q(t),
    D'(t) = k*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = beta*I(t)*S(t)/N(t)-alpha*S(t),
    E'(t) =  beta*I(t)*S(t)/N(t)-gamma*E(t),
    I'(t) = gamma*E(t)-delta*I(t),
    R'(t) = delta*I(t)-u1(t)*Q(t)-u2(t)*Q(t),
    Q'(t) = u1(t)*Q(t),
    D'(t) = u2(t)*Q(t),
    P'(t) = alpha*S(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

