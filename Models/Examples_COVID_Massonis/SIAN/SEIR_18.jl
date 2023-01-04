#Modeling the Control of COVID-19: Impact of Policy Interventions and Meteorological Factors

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*(I(t)+theta*A(t))-pp*S(t)+lambda*Q(t),
    Q'(t) = pp*S(t)-lambda*Q(t),
    E'(t) = beta*S(t)*(I(t)+theta*A(t))-sigma*E(t),
    A'(t) = sigma*(1-rho)*E(t)-e_a*A(t)-gamma_a*A(t),
    I'(t) = sigma*rho*E(t)-gamma_i*I(t)-d_i*I(t)-e_i*I(t),
    D'(t) = e_a*A(t)+e_i*I(t)-d_d*D(t)-gamma_d*D(t),
    R'(t) = gamma_a*A(t) + gamma_i*I(t) + gamma_d*D(t),
    y1(t) = D(t),
    y2(t) = Q(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#
using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -beta*S(t)*(I(t)+theta*A(t))-pp*S(t)+lambda*Q(t),
    Q'(t) = pp*S(t)-lambda*Q(t),
    E'(t) = beta*S(t)*(I(t)+theta*A(t))-sigma*E(t),
    A'(t) = sigma*(1-rho)*E(t)-e_a*A(t)-gamma_a*A(t),
    I'(t) = sigma*rho*E(t)-gamma_i*I(t)-d_i*I(t)-e_i*I(t),
    D'(t) = e_a*A(t)+e_i*I(t)-d_d*D(t)-gamma_d*D(t),
    R'(t) = gamma_a*A(t) + gamma_i*I(t) + gamma_d*D(t),
    y1(t) = D(t),
    y2(t) = Q(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

