# SEIR Transmission dynamics model of 2019 nCoV coronavirus with considering the 
# weak infectious ability and changes in latency duration


using SIAN, Logging
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

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#PRUEBAS
#todos SGI
using SIAN, Logging
ode = @ODEmodel(
    sigma'(t) = 0,
    q'(t) = 0,
    S'(t) = -(c*be+c+q(t)*(1-be))*S(t)*(I(t)+theta*E(t))+la*Sq(t),
    E'(t) = c*be*(1-q(t))*S(t)*(I(t)+theta*E(t))-sigma(t)*E(t),
    I'(t) = sigma(t)*E(t)-(delta_i+alpha+gamma_i)*I(t),
    Sq'(t) = c*q(t)*(1-be)*S(t)*(I(t)+theta*E(t))-la*Sq(t),
    Eq'(t) = c*be*q(t)*S(t)*(I(t)+theta*E(t))-Eq(t),#eliminé deltaQ
    H'(t) = delta_i*I(t)+delta_q*Eq(t)-(alpha+gamma_h)*H(t),
    R'(t) = gamma_i*I(t)+gamma_h*H(t),
    y1(t) = I(t)*vi,
    y2(t) = R(t)*vr,
    y3(t) = sigma(t),
    y4(t) = q(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#todo sigue igual
using SIAN, Logging
ode = @ODEmodel(
    sigma'(t) = 0,
    q'(t) = 0,
    S'(t) = -(c*be+c+q(t)*(1-be))*S(t)*(I(t)+theta*E(t))+la*Sq(t),
    E'(t) = c*be*(1-q(t))*S(t)*(I(t)+theta*E(t))-sigma(t)*E(t),
    I'(t) = sigma(t)*E(t)-(delta_i+alpha+gamma_i)*I(t),
    Sq'(t) = c*q(t)*(1-be)*S(t)*(I(t)+theta*E(t))-la*Sq(t),
    Eq'(t) = c*be*q(t)*S(t)*(I(t)+theta*E(t))-delta_q*Eq(t),
    H'(t) = delta_i*I(t)+delta_q*Eq(t)-(alpha+gamma_h)*H(t),
    R'(t) = gamma_i*I(t)+H(t),#ELIMINÉ GAMMAH
    y1(t) = I(t)*vi,
    y2(t) = R(t)*vr,
    y3(t) = sigma(t),
    y4(t) = q(t)
)


@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#
using SIAN, Logging
ode = @ODEmodel(
    sigma'(t) = 0,
    q'(t) = 0,
    S'(t) = -(c*be+c+q(t)*(1-be))*S(t)*(I(t)+theta*E(t))+la*Sq(t),
    E'(t) = c*be*(1-q(t))*S(t)*(I(t)+theta*E(t))-sigma(t)*E(t),
    I'(t) = sigma(t)*E(t)-(delta_i+alpha+gamma_i)*I(t),
    Sq'(t) = c*q(t)*(1-be)*S(t)*(I(t)+theta*E(t))-la*Sq(t),
    Eq'(t) = c*be*q(t)*S(t)*(I(t)+theta*E(t))-delta_q*Eq(t),
    H'(t) = delta_i*I(t)+Eq(t)-(alpha+gamma_h)*H(t),
    R'(t) = gamma_i*I(t)+gamma_h*H(t),
    y1(t) = I(t)*vi,
    y2(t) = R(t)*vr,
    y3(t) = sigma(t),
    y4(t) = q(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
