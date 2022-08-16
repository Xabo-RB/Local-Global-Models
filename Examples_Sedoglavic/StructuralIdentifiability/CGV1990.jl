using StructuralIdentifiability

ode = @ODEmodel(
    q1'(t) = k4*q3(t)-(k3+k7)*q1(t) + u(t),
    q3'(t) = k3*q1(t)-k4*q3(t) - k5*q3(t)*(R*V3-q35(t))+k6*q35(t)-k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))+k6*q36(t),
    q35'(t) = k5*q3(t)*(R*V3-q35(t))-k6*q35(t),
    q36'(t) = k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))-k6*q36(t),
    q7'(t) = k7*q1(t),
    y1(t) = q7(t)
)

@time println(assess_identifiability(ode))


#SIN u(t)
using StructuralIdentifiability

ode = @ODEmodel(
    q1'(t) = k4*q3(t)-(k3+k7)*q1(t),
    q3'(t) = k3*q1(t)-k4*q3(t) - k5*q3(t)*(R*V3-q35(t))+k6*q35(t)-k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))+k6*q36(t),
    q35'(t) = k5*q3(t)*(R*V3-q35(t))-k6*q35(t),
    q36'(t) = k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))-k6*q36(t),
    q7'(t) = k7*q1(t),
    y1(t) = q7(t)
)

@time println(assess_identifiability(ode))

#####################################################################################
using SIAN, Logging
ode = @ODEmodel(
    q1'(t) = k4*q3(t)-(k3+k7)*q1(t)+u(t),
    q3'(t) = k3*q1(t)-k4*q3(t) - k5*q3(t)*(R*V3-q35(t))+k6*q35(t)-k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))+k6*q36(t),
    q35'(t) = k5*q3(t)*(R*V3-q35(t))-k6*q35(t),
    q36'(t) = k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))-k6*q36(t),
    q7'(t) = k7*q1(t),
    y1(t) = q7(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


using SIAN, Logging
ode = @ODEmodel(
    q1'(t) = k4*q3(t)-(k3+k7)*q1(t),
    q3'(t) = k3*q1(t)-k4*q3(t) - k5*q3(t)*(R*V3-q35(t))+k6*q35(t)-k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))+k6*q36(t),
    q35'(t) = k5*q3(t)*(R*V3-q35(t))-k6*q35(t),
    q36'(t) = k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))-k6*q36(t),
    q7'(t) = k7*q1(t),
    y1(t) = q7(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))