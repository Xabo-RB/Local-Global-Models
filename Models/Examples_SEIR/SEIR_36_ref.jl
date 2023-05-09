# https://github.com/ryansmcgee/seirsplus/blob/master/docs/SEIRSplus_Model.pdf

using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    nu'(t) = 0,
    q'(t) = 0,
    S'(t) = -beta*S(t)*I(t)/N(t)-q(t)*beta_d*S(t)*Di(t)/N(t)+nu(t)*N(t)-mu_0*S(t),
    E'(t) = beta*S(t)*I(t)/N(t)+q(t)*beta_d*S(t)*Di(t)/N(t)-s*E(t)-phi_e*E(t)-mu_0*E(t),
    I'(t) = s*E(t)-gamma*I(t)-mu_i*I(t)-phi*I(t)-mu_0*I(t),
    De'(t) = phi_e*E(t)-s_d*De(t)-mu_0*De(t),
    Di'(t) = phi*I(t)+s_d*De(t)-gamma_d*Di(t)-mu_d*Di(t)-mu_0*Di(t),
    R'(t) = gamma*I(t)+gamma_d*Di(t)-mu_0*R(t),
    F'(t) = mu_i*I(t)+mu_d*Di(t),
    y1(t) = De(t),
    y2(t) = Di(t),
    y5(t) = F(t),
    y3(t) = N(t),
    y4(t) = nu(t),
    y5(t) = q(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
