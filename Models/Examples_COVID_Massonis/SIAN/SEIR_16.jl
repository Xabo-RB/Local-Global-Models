#Why is it difficult to accurately predict the COVID-19 epidemic?


using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*I(t)*S(t),
    E'(t) = beta*I(t)*S(t)-epsilon*E(t),
    I'(t) = epsilon*E(t)-(rho+mu)*I(t),
    R'(t) = rho*I(t)-d*R(t),
    y1(t) = mu*I(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
