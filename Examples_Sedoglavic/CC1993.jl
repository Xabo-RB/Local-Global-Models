using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -(kp+(F01/V1)/g+k21)*x1(t) + k12*x2(t),
    x2'(t) = k21*x1(t)-(k02+x3(t)+k12)*x2(t),
    x3'(t) = -kb*x3(t)+ka*u(t),
    y1(t) = x1(t)/V1
)

@time println(assess_identifiability(ode))
