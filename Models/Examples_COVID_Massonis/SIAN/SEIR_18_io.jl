using StructuralIdentifiability


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

#@time println(assess_identifiability(ode))

println(find_ioequations(ode))
