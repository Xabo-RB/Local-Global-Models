using StructuralIdentifiability

ode = @ODEmodel(
    q1'(t) = k4*q3(t)-(k3+k7)*q1(t) + u(t),
    q3'(t) = k3*q1(t)-k4*q3(t) - k5*q3(t)*(R*V3-q35(t))+k6*q35(t)-k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))+k6*q36(t),
    q35'(t) = k5*q3(t)*(R*V3-q35(t))-k6*q35(t),
    q36'(t) = k5*q3(t)*(5*V36/V3)*(S*V36-q36(t))-k6*q36(t),
    q7'(t) = k7*q1(t),
    y1(t) = q7(t)
)

println(find_ioequations(ode))
