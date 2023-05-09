# Construction of Compartmental Models for COVID-19 with Quarantine, Lockdown and Vaccine Interventions


using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    mu'(t) = 0,
    S'(t) = mu(t)*N(t)-beta_1*S(t)*I(t)-(gamma+eta)*S(t)+delta*L(t)+xi*E(t),
    L'(t) = eta*S(t)-(gamma+delta)*L(t),
    E'(t) = beta_1*S(t)*I(t)-(gamma+theta_2+epsilon+xi)*E(t),
    I'(t) = epsilon*E(t)-(gamma+theta_1+alpha_1)*I(t),
    Q'(t) = theta_1*I(t)+theta_2*E(t)-(gamma+alpha_2)*Q(t),
    R'(t) = alpha_1*I(t)+alpha_2*Q(t)-gamma*R(t),
    y1(t) = L(t),
    y2(t) = Q(t),
    y3(t) = N(t),
    y4(t) = mu(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

