# Modelling the transmission dynamics of COVID-19
# in six high burden countries


using SIAN, Logging

ode = @ODEmodel(
    A'(t) = 0,
    S'(t) = A(t)-beta*S(t)*I(t)/(1+alpha*I(t))-mu*S(t),
    L'(t) = beta*S(t)*I(t)/(1+alpha*I(t))-(W+mu)*L(t),
    I'(t) = W*L(t)-(gamma+mu)*I(t),
    R'(t) = gamma*I(t)-mu*R(t),
    y1(t) = W*L(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
