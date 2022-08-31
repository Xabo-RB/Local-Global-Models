#IDENTIFIABILITY OF INFECTION MODEL PARAMETERS EARLY IN AN EPIDEMIC. TIMOTHY SAUER 2021

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta* (U(t)+I(t)) *(S(t)/N),
    E'(t) =  beta*(U(t)+I(t))*(S(t)/N) - E(t)/Z,
    U'(t) = (1-alpha)*E(t)/Z - U(t)/D,
    U'(t) = alpha*E(t)/Z - I(t)/D,
    R'(t) = U(t)/D + I(t)/D,
    y1(t) = I(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
