# Research about the optimal strategies for prevention and control 
# of varicella outbreak in a school in a central city of China: based on an SEIR dynamic model

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t),
    E'(t) =  beta*S(t)*I(t)-v*E(t),
    I'(t) = v*E(t)-psi*I(t)-(1-psi)*gamma*I(t),
    R'(t) = gamma*Q(t)+(1-psi)*gamma*I(t),
    Q'(t) = -gamma*Q(t)+psi*I(t),
    y1(t) = Q(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#PRUEBAS
#
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = beta*S(t)*I(t), #quité el signo
    E'(t) = beta*S(t)*I(t)-v*E(t),
    I'(t) = v*E(t)-psi*I(t)-(1-psi)*gamma*I(t),
    R'(t) = gamma*Q(t)+(1-psi)*gamma*I(t),
    Q'(t) = -gamma*Q(t)+psi*I(t),
    y1(t) = Q(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = beta*S(t)*I(t), #quité el signo
    E'(t) = beta*S(t)*I(t)-v*E(t),
    I'(t) = -v*E(t)-psi*I(t)-(1-psi)*gamma*I(t),#añadí el signo
    R'(t) = gamma*Q(t)+(1-psi)*gamma*I(t),
    Q'(t) = -gamma*Q(t)+psi*I(t),
    y1(t) = Q(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t), 
    E'(t) = beta*S(t)*I(t)-v*E(t),
    I'(t) = v + E(t)-psi*I(t)-(1-psi)*gamma*I(t), #deshice simetría
    R'(t) = gamma*Q(t)+(1-psi)*gamma*I(t),
    Q'(t) = -gamma*Q(t)+psi*I(t),
    y1(t) = Q(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
