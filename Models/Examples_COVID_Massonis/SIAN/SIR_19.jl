# Gaeta, G. (2020). A simple SIR model with a large set of asymptomatic infectives. arXiv preprint arXiv:2003.08720.

using SIAN, Logging
ode = @ODEmodel(
    S'(t) = -alpha*(I(t)+J(t))*S(t),
    I'(t) = alpha*xi*S(t)*(I(t)+J(t))-beta*I(t),
    J'(t) = alpha*(1-xi)*S(t)*(I(t)+J(t))-eta*J(t),
    R'(t) = beta*I(t),
    U'(t) = eta*J(t),
    y1(t) = R(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
