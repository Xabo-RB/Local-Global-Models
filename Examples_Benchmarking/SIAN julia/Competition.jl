
using Logging, SIAN


ode = @ODEmodel(
    x1'(t) = x1*(a1-b11*log(x1(t))-b12*log(x2(t))),
    x2'(t) = x2(t)*(a2-b21*log(x1(t))-b22*log(x2(t))),
    y1(t) = x1
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
