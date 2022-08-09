using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    X'(t) = mumax*S(t)*X(t)/(S(t)+Ks)-d*X(t),
    S'(t) = -mumax*S(t)*X(t)/(S(t)+Ks)/Y + d*(Sin - S(t)) ,
    y1(t) = S(t)
)

@time println(assess_identifiability(ode))

