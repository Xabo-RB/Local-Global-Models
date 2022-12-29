#Panteleev, M. A., Balandina, A. N., Lipets, E. N., Ovanesov, M. V., & Ataullakhanov, F. I. (2010). Task-oriented modular decomposition of biological networks: trigger mechanism in blood coagulation.
# Biophysical journal, 98(9), 1751-1761.

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta* (U(t)+I(t)) *(S(t)/N),
    E'(t) =  beta*(U(t)+I(t))*(S(t)/N) - E(t)/Z,
    U'(t) = (1-alpha)*E(t)/Z - U(t)/D,
    I'(t) = alpha*E(t)/Z - I(t)/D,
    R'(t) = U(t)/D + I(t)/D,
    y1(t) = I(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
