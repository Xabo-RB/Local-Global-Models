using StructuralIdentifiability

ode = @ODEmodel(
    xA'(t) = -k1 * xA(t),
    xB'(t) = k1 * xA(t) - k2 * xB(t),
    xC'(t) = k2 * xB(t),
    eA'(t) = 0,
    eC'(t) = 0,
    y1(t) = xC(t),
    y2(t) = eA(t) * xA(t) + eB * xB(t) + eC(t) * xC(t),
    y3(t) = eA(t),
    y4(t) = eC(t)
)


println(find_ioequations(ode))
