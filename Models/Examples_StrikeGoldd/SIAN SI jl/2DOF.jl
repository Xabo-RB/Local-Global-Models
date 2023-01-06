using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = dx1(t),
    x2'(t) = dx2(t),
    dx1(t) = 1/m1*( -(k1+dk1*x1(t))*x1(t)+k2*(x2(t)-x1(t))-c1*dx1(t)+c2*(dx2(t)-dx1(t))+F1 ),
    dx2(t) = 1/m2*( k2*(x1(t)-x2(t))+c2*(dx1(t)-dx2(t))+F2),
    y1(t) = x1(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#PRUEBAS
#No cambia nada
using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = dx1(t),
    x2'(t) = dx2(t),
    dx1(t) = 1/m1*( -(k1+dk1*x1(t))*x1(t)+k2*(x2(t)-x1(t))-c1*dx1(t)+(dx2(t)-dx1(t))+F1 ),#elimino c2
    dx2(t) = 1/m2*( k2*(x1(t)-x2(t))+c2*(dx1(t)-dx2(t))+F2),
    y1(t) = x1(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#Se pierde identificabilidad x1 baja a SLI y c2 baja a NI
using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = dx1(t),
    x2'(t) = dx2(t),
    dx1(t) = 1/m1*( -(k1+dk1*x1(t))*x1(t)+k2*(x2(t)-x1(t))-c1*dx1(t)+c2*(dx2(t)-dx1(t))+F1 ),
    dx2(t) = 1/m2*( k2*(x1(t)-x2(t))+(dx1(t)-dx2(t))+F2),
    y1(t) = x1(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
