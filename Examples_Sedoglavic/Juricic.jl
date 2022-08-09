using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -1/(Ci*Ri)*x1(t) + 1/(Ci*Ri)*x(t),
    x2'(t) = -1/Cw*(1/Rw1 +1/Rw2)*x2(t) + 1/(Cw*Rw1)*x3(t) + 1/(Cw*Rw2)*u1(t),
    q35'(t) = 1/(Ca*Ri)*x1(t) + 1/(Ca*Rw1)*x2(t) - 1/Ca*(1/Ri + 1/Rw1 + zetaD/RD + 1/Rn + zetaW/RW +1/Rout)*x3(t) + 1/Ca*(zetaD/RD+1/Rn)*u1(t) + 1/Ca*(zetaW/RW +1/Rout)*u2(t) + 1/Ca*u3(t),
    y1(t) = x3(t)
)

@time println(assess_identifiability(ode))