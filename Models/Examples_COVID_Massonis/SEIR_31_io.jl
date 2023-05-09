using StructuralIdentifiability


ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -alp*S(t)*(I(t)+E(t)+A(t))/N(t)-sig*N(t),
    E'(t) = alp*S(t)*(I(t)+E(t)+A(t))/N(t)-beta_1*E(t),
    I'(t) = beta_1*hh*E(t)+beta_2*r*A(t)-phi*q*I(t)-gamma*(1-q)*I(t),
    A'(t) = beta_1*(1-hh)*E(t)-beta_2*r*A(t)-gamma*(1-r)*A(t),
    J'(t) = phi*q*I(t)-gamma*J(t),
    R'(t) = gamma*(1-q)*I(t)+gamma*(1-r)*A(t)+gamma*J(t),
    y1(t) = I(t),
    y2(t) = J(t),
    y3(t) = R(t),
    y4(t) = N(t)
)

#@time println(assess_identifiability(ode))

println(find_ioequations(ode))
