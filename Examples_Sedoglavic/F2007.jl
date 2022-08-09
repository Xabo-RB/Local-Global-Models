using StructuralIdentifiability

ode = @ODEmodel(
    A1'(t) = -(k10+k12+k13)*A1(t)+k21*A2(t)+k31*A3(t)+u(t),
    A2'(t) = k12*A1(t)-k21*A2(t),
    A3'(t) = k13*A1(t)-k31*A3(t),
    y1(t) = (1+B1*u2(t)+B2+B3/(K+A1(t)/V1))*A1(t)/V1
)

@time println(assess_identifiability(ode))