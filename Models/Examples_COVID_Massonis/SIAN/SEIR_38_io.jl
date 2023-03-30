using StructuralIdentifiability


ode = @ODEmodel(
    sigma'(t) = 0,
    q'(t) = 0,
    S'(t) = -(c*be+c+q(t)*(1-be))*S(t)*(I(t)+theta*E(t))+la*Sq(t),
    E'(t) = c*be*(1-q(t))*S(t)*(I(t)+theta*E(t))-sigma(t)*E(t),
    I'(t) = sigma(t)*E(t)-(delta_i+alpha+gamma_i)*I(t),
    Sq'(t) = c*q(t)*(1-be)*S(t)*(I(t)+theta*E(t))-la*Sq(t),
    Eq'(t) = c*be*q(t)*S(t)*(I(t)+theta*E(t))-delta_q*Eq(t),
    H'(t) = delta_i*I(t)+delta_q*Eq(t)-(alpha+gamma_h)*H(t),
    R'(t) = gamma_i*I(t)+gamma_h*H(t),
    y1(t) = I(t)*vi,
    y2(t) = R(t)*vr,
    y3(t) = sigma(t),
    y4(t) = q(t)
)

#@time println(assess_identifiability(ode))

println(find_ioequations(ode))
