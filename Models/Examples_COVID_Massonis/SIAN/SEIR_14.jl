# A modified SEIR model to predict the COVID-19 outbreak in Spain
# and Italy: simulating control scenarios and multi-scale epidemics

using SIAN, Logging

ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = mu*N(t) - alpha*S(t) - beta*I(t)*S(t) - mu*S(t),
    E'(t) = beta*I(t)*S(t)*N(t)-mu*E(t)-gamma*E(t),
    I'(t) = gamma*E(t)-delta*I(t)-mu*I(t)-mu*S(t),
    Q'(t) = delta*I(t)-(1-lambda_0*exp(-lambda_1*u2(t)))*Q(t)-(k_0*exp(-k_1*u1(t)))*Q(t)-mu*Q(t),
    R'(t) = (1-lambda_0*exp(-lambda_1*u2(t)))*Q(t)-mu*R(t),
    D'(t) = (k_0*exp(-k_1*u1(t)))*Q(t),
    C'(t) = alpha*S(t)-mu*C(t)-tau*C(t),
    y1(t) = Q(t),
    y2(t) = D(t),
    y3(t) = C(t),
    y4(t) = N(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
