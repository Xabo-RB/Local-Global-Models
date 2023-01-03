# To mask or not to mask: Modeling the potential for face mask
# use by the general public to curtail the COVID-19 pandemic

using SIAN, Logging

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

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#PRUEBAS
#no cambia nada
using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*(I(t)+eta*A(t))*S(t)/N(t),
    E'(t) = -beta*(I(t)+eta*A(t))*S(t)/N(t)-sigma*E(t), #AÑADÍ SIGNO
    I'(t) = alpha*sigma*E(t)-phi*I(t)-gamma_i*I(t),
    A'(t) = (1-alpha)*sigma*E(t)-gamma_0*A(t),
    H'(t) = phi*I(t)-delta*H(t)-gamma_h*H(t),
    R'(t) = gamma_i*I(t)+gamma_0*A(t)+gamma_h*H(t),
    D'(t) = delta*H(t),
    y1(t) = I(t),
    y2(t) = H(t),
    y3(t) = D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#ningun SLI, beta y sigma van a SGI y gammai a NI
using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*(I(t)+eta*A(t))*S(t)/N(t),
    E'(t) = beta*(I(t)+eta*A(t))*S(t)/N(t)-sigma*E(t), 
    I'(t) = alpha*sigma*E(t)-phi*I(t),#-gamma_i*I(t), ELIMINADO
    A'(t) = (1-alpha)*sigma*E(t)-gamma_0*A(t),
    H'(t) = phi*I(t)-delta*H(t)-gamma_h*H(t),
    R'(t) = gamma_i*I(t)+gamma_0*A(t)+gamma_h*H(t),
    D'(t) = delta*H(t),
    y1(t) = I(t),
    y2(t) = H(t),
    y3(t) = D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#NO INFLUYE
using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*(I(t)+eta*A(t))*S(t)/N(t),
    E'(t) = beta*(I(t)+eta*A(t))*S(t)/N(t)-sigma*E(t), 
    I'(t) = -phi*I(t)-gamma_i*I(t)+alpha*E(t),#QUITÉ SIGMA, SIGMA = 1
    A'(t) = (1-alpha)*sigma*E(t)-gamma_0*A(t),
    H'(t) = phi*I(t)-delta*H(t)-gamma_h*H(t),
    R'(t) = gamma_i*I(t)+gamma_0*A(t)+gamma_h*H(t),
    D'(t) = delta*H(t),
    y1(t) = I(t),
    y2(t) = H(t),
    y3(t) = D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#Se eliminan los SLI, beta y gammai a SGI y sigma a NI
using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*(I(t)+eta*A(t))*S(t)/N(t),
    E'(t) = beta*(I(t)+eta*A(t))*S(t)/N(t)-E(t),#QUITÉ SIGMA, SIGMA = 1
    I'(t) = alpha*E(t)-phi*I(t)-gamma_i*I(t),#QUITÉ SIGMA, SIGMA = 1
    A'(t) = (1-alpha)*sigma*E(t)-gamma_0*A(t),
    H'(t) = phi*I(t)-delta*H(t)-gamma_h*H(t),
    R'(t) = gamma_i*I(t)+gamma_0*A(t)+gamma_h*H(t),
    D'(t) = delta*H(t),
    y1(t) = I(t),
    y2(t) = H(t),
    y3(t) = D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#beta gammai pasan a SGI y sigma a NI
using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -beta*(I(t)+eta*A(t))*S(t)/N(t),
    E'(t) = (I(t)+eta*A(t))*S(t)/N(t)-sigma*E(t), #eliminé beta
    I'(t) = alpha*sigma*E(t)-phi*I(t)-gamma_i*I(t),
    A'(t) = (1-alpha)*sigma*E(t)-gamma_0*A(t),
    H'(t) = phi*I(t)-delta*H(t)-gamma_h*H(t),
    R'(t) = gamma_i*I(t)+gamma_0*A(t)+gamma_h*H(t),
    D'(t) = delta*H(t),
    y1(t) = I(t),
    y2(t) = H(t),
    y3(t) = D(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
