#IDENTIFIABILITY OF INFECTION MODEL PARAMETERS EARLY IN AN EPIDEMIC. TIMOTHY SAUER 2021
#The same as SEIR from SIANjl

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*I(t)*(S(t)/N(t)),
    E'(t) =  beta*I(t)*(S(t)/N(t)) - alpha*E(t),
    I'(t) = alpha*E(t) - lambda*I(t),
    R'(t) = lambda*I(t),
    y1(t) = I(t),
    y2(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
