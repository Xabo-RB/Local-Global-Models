# To mask or not to mask: Modeling the potential for face mask
# use by the general public to curtail the COVID-19 pandemic

using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*(I(t)+eta*A(t))*S(t)/N,
    E'(t) = beta*(I(t)+eta*A(t))*S(t)/N-sigma*E(t),
    I'(t) = alpha*sigma*E(t)-phi*I(t)-gamma_i*I(t),
    A'(t) = gamma*(1-pp)*E(t)-mu2*A(t),
    R'(t) = mu1*I(t)+mu2*A(t),
    y1(t) = I(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

