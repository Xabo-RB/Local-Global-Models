#Roda, W. C., Varughese, M. B., Han, D., & Li, M. Y. (2020). Why is it difficult to accurately predict the COVID-19 epidemic?.
# Infectious disease modelling, 5, 271-281.


using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*S(t)*I(t)/N(t)-pp*S(t)+q*C(t),
    I'(t) = beta*S(t)*I(t)/N(t)-(r+mu)*I(t),
    R'(t) = r*I(t),
    C'(t) = pp*S(t)-q*C(t),
    D'(t) = mu*I(t),
    y1(t) = N(t),
    y2(t) = D(t),
    y3(t) = C(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
