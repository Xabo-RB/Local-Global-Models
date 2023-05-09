using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -(k01+k21)*x1(t) + k12*x2(t) +k13*x1(t)*x3(t),
    x2'(t) = k21*x1(t) -(k02+k12+k32+k42)*x2(t) - k23*x3(t)*x2(t) + k24*x4(t),
    x3'(t) = k32*x2(t) - (k03+k13*x1(t)+k23*x2(t)-k43*x4(t))*x3(t)+k34*x4(t)+u(t),
    x4'(t) = k42*x2(t)+k43*x3(t)*x4(t)-(k24+k34)*x4(t),
    y1(t) = x1(t)
)

@time println(assess_identifiability(ode))