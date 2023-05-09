# Gaeta, G. (2020). A simple SIR model with a large set of asymptomatic infectives. arXiv preprint arXiv:2003.08720.

using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -alpha*(I(t)+J(t))*S(t),
    I'(t) = alpha*xi*S(t)*(I(t)+J(t))-beta*I(t),
    J'(t) = alpha*(1-xi)*S(t)*(I(t)+J(t))-eta*J(t),
    R'(t) = beta*I(t),
    U'(t) = eta*J(t),
    y1(t) = R(t)
)
#@time println(assess_identifiability(ode))

println(find_ioequations(ode))
