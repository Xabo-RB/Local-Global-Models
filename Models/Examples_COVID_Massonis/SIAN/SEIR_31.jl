# Mathematical model of transmission dynamics with mitigation and health measures for SARS-CoV-2 infection in European countries

using SIAN, Logging
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
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#PRUEBAS
#todos SGI
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -alp*S(t)*(I(t)+E(t)+A(t))/N(t)-sig*N(t),
    E'(t) = alp*S(t)*(I(t)+E(t)+A(t))/N(t)-E(t),#ELIMINO BETA1
    I'(t) = beta_1*hh*E(t)+beta_2*r*A(t)-phi*q*I(t)-gamma*(1-q)*I(t),
    A'(t) = beta_1*(1-hh)*E(t)-beta_2*r*A(t)-gamma*(1-r)*A(t),
    J'(t) = phi*q*I(t)-gamma*J(t),
    R'(t) = gamma*(1-q)*I(t)+gamma*(1-r)*A(t)+gamma*J(t),
    y1(t) = I(t),
    y2(t) = J(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#todos SGI
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -alp*S(t)*(I(t)+E(t)+A(t))/N(t)-sig*N(t),
    E'(t) = alp*S(t)*(I(t)+E(t)+A(t))/N(t)-beta_1*E(t),
    I'(t) = hh*E(t)+beta_2*r*A(t)-phi*q*I(t)-gamma*(1-q)*I(t),#ELIMINO BETA1
    A'(t) = beta_1*(1-hh)*E(t)-beta_2*r*A(t)-gamma*(1-r)*A(t),
    J'(t) = phi*q*I(t)-gamma*J(t),
    R'(t) = gamma*(1-q)*I(t)+gamma*(1-r)*A(t)+gamma*J(t),
    y1(t) = I(t),
    y2(t) = J(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#todo sgi
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -alp*S(t)*(I(t)+E(t)+A(t))/N(t)-sig*N(t),
    E'(t) = alp*S(t)*(I(t)+E(t)+A(t))/N(t)-beta_1*E(t),
    I'(t) = beta_1*hh*E(t)+beta_2*r*A(t)-phi*q*I(t)-gamma*(1-q)*I(t),
    A'(t) = (1-hh)*E(t)-beta_2*r*A(t)-gamma*(1-r)*A(t),
    J'(t) = phi*q*I(t)-gamma*J(t),
    R'(t) = gamma*(1-q)*I(t)+gamma*(1-r)*A(t)+gamma*J(t),
    y1(t) = I(t),
    y2(t) = J(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#todos SGI
using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -alp*S(t)*(I(t)+E(t)+A(t))/N(t)-sig*N(t),
    E'(t) = alp*S(t)*(I(t)+E(t)+A(t))/N(t)-beta_1*E(t),
    I'(t) = beta_1*hh*E(t)+r*A(t)-phi*q*I(t)-gamma*(1-q)*I(t),#quito beta2
    A'(t) = beta_1*(1-hh)*E(t)-beta_2*r*A(t)-gamma*(1-r)*A(t),
    J'(t) = phi*q*I(t)-gamma*J(t),
    R'(t) = gamma*(1-q)*I(t)+gamma*(1-r)*A(t)+gamma*J(t),
    y1(t) = I(t),
    y2(t) = J(t),
    y3(t) = R(t),
    y4(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
