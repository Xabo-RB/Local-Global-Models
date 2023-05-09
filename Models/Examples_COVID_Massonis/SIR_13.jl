# http://rocs.hu-berlin.de/corona/docs/forecast/model/

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -alpha*I(t)*S(t)-k0*S(t),
    I'(t) = alpha*I(t)*S(t)-beta*I(t)-k*I(t)-k0*S(t),
    X'(t) = (k+k0)*I(t),
    R'(t) = beta*I(t)+k0*S(t),
    y1(t) = X(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
