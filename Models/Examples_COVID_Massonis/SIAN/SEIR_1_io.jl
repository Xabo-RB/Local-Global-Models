# Research about the optimal strategies for prevention and control 
# of varicella outbreak in a school in a central city of China: based on an SEIR dynamic model

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)


ode = @ODEmodel(
    S'(t) = -beta*S(t)*I(t),
    E'(t) =  beta*S(t)*I(t)-v*E(t),
    I'(t) = v*E(t)-psi*I(t)-(1-psi)*gamma*I(t),
    R'(t) = gamma*Q(t)+(1-psi)*gamma*I(t),
    Q'(t) = -gamma*Q(t)+psi*I(t),
    y1(t) = Q(t)
)
@time println(assess_identifiability(ode))

find_ioequations(ode, [var_change_policy=:default])

