using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = dx1(t),
    x2'(t) = dx2(t),
    dx1(t) = 1/m1*( -(k1+dk1*x1(t))*x1(t)+k2*(x2(t)-x1(t))-c1*dx1(t)+c2*(dx2(t)-dx1(t))+F1 ),
    dx2(t) = 1/m2*( k2*(x1(t)-x2(t))+c2*(dx1(t)-dx2(t))+F2),
    y1(t) = x1(t)
)

println(find_ioequations(ode))
