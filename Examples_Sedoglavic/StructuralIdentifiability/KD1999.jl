using StructuralIdentifiability

ode = @ODEmodel(
    Ca'(t) = u1(t)*(Ca0-Ca(t))/V - k0*Arr*Ca(t),
    Cb'(t) = -u1(t)*Cb(t)/V + k0*Arr(t)*Ca(t),
    T'(t) = u1(t)*(Ta-T(t))/V - (k0*Arr(t)*Ca(t)*DH + UA*(Tj(t)-T(t))/V)/(ro*cp),
    Tj'(t) = u2(t)*(Th-Tj(t))/Vh - UA/(roh*cph)*(Tj(t)-T(t))/Vh,
    Arr'(t) = E*Arr(t)/(R*T(t)^2)*(u1(t)*(Ta-T(t))/V - (k0*Arr(t)*Ca(t)*DH + UA*(Tj(t)-T(t))/V)/(ro*cp)),
    y1(t) = Cb(t),
    y2(t) = T(t)
)

@time println(assess_identifiability(ode))