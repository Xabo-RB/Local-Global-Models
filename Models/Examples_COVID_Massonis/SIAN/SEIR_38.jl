# SEIR Transmission dynamics model of 2019 nCoV coronavirus with considering the 
# weak infectious ability and changes in latency duration


using SIAN, Logging
ode = @ODEmodel(
    sigma'(t) = 0,
    q'(t) = 0,
    S'(t) = -(c*be+c+q*(1-be))*S(t)*(I(t)+theta*E(t))+la*Sq(t),
    E'(t) = c*be*(1-q)*S(t)*(I(t)+theta*E(t))-sigma*E(t),
    I'(t) = sigma*E(t)-(delta_i+alpha+gamma_i)*I(t),
    Sq'(t) = c*q*(1-be)*S(t)*(I(t)+theta*E(t))-la*Sq(t),
    Eq'(t) = c*be*q*S(t)*(I(t)+theta*E(t))-delta_q*Eq(t),
    H'(t) = delta_i*I(t)+delta_q*Eq(t)-(alpha+gamma_h)*H(t),
    R'(t) = gamma_i*I(t)+gamma_h*H(t),
    y1(t) = I(t)*vi,
    y2(t) = R(t)*vr,
    y3(t) = sigma(t),
    y4(t) = q(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
