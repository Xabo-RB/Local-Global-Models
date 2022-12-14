#Kim, M. Y., & Milner, F. A. (1995). A mathematical model of epidemics with screening and variable infectivity.
#Mathematical and computer modelling, 21(7), 29-42.


using SIAN, Logging
ode = @ODEmodel(
    A'(t) = 0,
    S'(t) = A(t)-mu+S(t)-c*phi*I(t)*S(t)/(I(t)+S(t)),
    I'(t) = -mu*I(t)+c*phi*I(t)*S(t)/(I(t)+S(t))-gamma*I(t)-I(t)*u1(t)/(S(t)+I(t)),
    R'(t) = -mu*R(t)+gamma*I(t)+I(t)*u1(t)/(S(t)+I(t)),
    y1(t) = K*I(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
