
#Primera reparametrización
using Logging, SIAN


ode = @ODEmodel(
    s'(t) = mu + g*r(t) - mu*s(t) - b0*i(t)*s(t) - b0*i(t)*s(t)*x1(t),
    x2'(t) = M*x1(t),
    i'(t) = -i(t)*(mu + nu - b0*s(t) - b0*s(t)*x1(t)),
    r'(t) = i(t)*nu - g*r(t) - mu*r(t),
    x1'(t) = -M*x2(t),
    y1(t) = i(t),
    y2(t) = r(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1))
