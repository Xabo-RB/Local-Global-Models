# Wenjie Zheng (2020) Total Variation Regularization for Compartmental Epidemic Models with Time-varying Dynamic

using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*I(t)*S(t)/N(t),
    I'(t) = beta*I(t)*S(t)/N(t)-gamma*I(t),
    R'(t) = gamma*I(t),
    y1(t) = I(t)*K,
    y2(t) = N(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

