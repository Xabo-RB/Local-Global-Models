# Mathematical model of transmission dynamics with mitigation and health measures for SARS-CoV-2 infection in European countries

using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -alp*S(t)*(I(t)+E(t)+A(t))/N(t)-sig*N(t),
    E'(t) = alp*S(t)*(I(t)+E(t)+A(t))/N(t)-beta_1*E(t),
    I'(t) = beta_1*hh*E(t)+beta_2*r*A(t)-phi*q*I(t)-gamma*(1-q)*I(t),
    A'(t) = beta_1*(1-hh)*E(t)-beta_2*r*A(t)-gamma*(1-r)*A(t),
    J'(t) = phi*q*I(t)-gamma*J(t),
    R'(t) = gamma*(1-q)*I(t)+gamma*(1-r)*A(t)+gamma*J(t),
    y1(t) = D(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
