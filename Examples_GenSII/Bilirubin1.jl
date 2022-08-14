using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -k31*x1(t)+k13*x3(t) + u(t),
    x2'(t) = -k42*x2(t)+k24*x4(t),
    x3'(t) = k31*x1(t)-(k03+k13+k43)*x3(t),
    x4'(t) = k42*x2(t)+k43*x3(t)-(k04+k24)*x4(t),
    y1(t) = x1(t),
    y2(t) = x2(t)
)

@time println(assess_identifiability(ode))