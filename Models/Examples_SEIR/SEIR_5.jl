# Mathematical Modeling of Epidemic Diseases; A Case Study of the COVID-19 Coronavirus

using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*(I(t)+A(t))*S(t),
    E'(t) = beta*(I(t)+A(t)*S(t)-gamma*E(t)),
    I'(t) = gamma*pp*E(t)-mu1*I(t),
    A'(t) = gamma*(1-pp)*E(t)-mu2*A(t),
    R'(t) = mu1*I(t)+mu2*A(t),
    y1(t) = I(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

