# Mathematical Modeling of Epidemic Diseases; A Case Study of the COVID-19 Coronavirus

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -alpha_e*S(t)*E(t)-alpha_i*S(t)*I(t),
    E'(t) =  alpha_e*S(t)*E(t)+alpha_i*S(t)*I(t)-k*E(t)-rho*E(t),
    I'(t) = k*E(t)-beta*I(t)-mu*I(t),
    R'(t) = beta*I(t)+rho*E(t),
    P'(t) = mu*I(t),
    y1(t) = P(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

