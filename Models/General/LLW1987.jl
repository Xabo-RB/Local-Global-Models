using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2*u(t),
    x2'(t) = - p3*x2(t) + p4*u(t),
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+p2*x2(t))*u(t),
    y1(t) = x3(t)
)

@time println(assess_identifiability(ode))


#SIN u(t)
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2,
    x2'(t) = - p3*x2(t) + p4,
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+p2*x2(t)),
    y1(t) = x3(t)
)

@time println(assess_identifiability(ode))

#CON otro u(t)
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2 + u(t),
    x2'(t) = - p3*x2(t) + p4 + u(t),
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+p2*x2(t)) +u(t),
    y1(t) = x3(t)
)

@time println(assess_identifiability(ode))


#CON otro u(t)
using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2 + u(t),
    x2'(t) = - p3*x2(t) + p4,
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+p2*x2(t)),
    y1(t) = x3(t)
)

@time println(assess_identifiability(ode))

#p1 y p3 a SGI
using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = - x1(t) + p2*u(t),#elimino p1
    x2'(t) = - p3*x2(t) + p4*u(t),
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+p2*x2(t))*u(t),
    y1(t) = x3(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#normal con SIAN 
using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2*u(t),
    x2'(t) = - p3*x2(t) + p4*u(t),
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+p2*x2(t))*u(t),
    y1(t) = x3(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#p1 y p3 SGI
using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2*u(t),
    x2'(t) = x2(t) + p4*u(t),#elimino p3
    x3'(t) = - (p1+p3)*x3(t) + (p4*x1(t)+p2*x2(t))*u(t),
    y1(t) = x3(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#p1 y p3 SGI
using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2*u(t),
    x2'(t) = -p3*x2(t) + p4*u(t),
    x3'(t) = - (p1)*x3(t) + (p4*x1(t)+p2*x2(t))*u(t),#elimino p3
    y1(t) = x3(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#p1 y p3 SGI
using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = - p1*x1(t) + p2*u(t),
    x2'(t) = -p3*x2(t) + p4*u(t),
    x3'(t) = - (p3)*x3(t) + (p4*x1(t)+p2*x2(t))*u(t),#elimino p1
    y1(t) = x3(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
