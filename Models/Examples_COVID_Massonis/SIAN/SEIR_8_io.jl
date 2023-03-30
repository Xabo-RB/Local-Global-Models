using StructuralIdentifiability


ode = @ODEmodel(
    S'(t) = -beta*(I(t)+eta*A(t))*S(t)/N(t),
    E'(t) = beta*(I(t)+eta*A(t))*S(t)/N(t)-sigma*E(t),
    I'(t) = alpha*sigma*E(t)-phi*I(t)-gamma_i*I(t),
    A'(t) = (1-alpha)*sigma*E(t)-gamma_0*A(t),
    H'(t) = phi*I(t)-delta*H(t)-gamma_h*H(t),
    R'(t) = gamma_i*I(t)+gamma_0*A(t)+gamma_h*H(t),
    D'(t) = delta*H(t),
    y1(t) = I(t),
    y2(t) = H(t),
    y3(t) = D(t)
)

#@time println(assess_identifiability(ode))

println(find_ioequations(ode))
