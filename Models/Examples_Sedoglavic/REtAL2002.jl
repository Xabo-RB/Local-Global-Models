using StructuralIdentifiability

ode = @ODEmodel(
    Qm'(t) = -kabs*Qm(t)+u(t),
    Qs'(t) = -((Ki+(Ks-kreabs)*(-(CCr0*(E50delta+Edelta(t))-Edelta(t)*CCrmax)/(-1+eta*(SC(t)*cosphi+CS(t)*sinphi))/(E50delta+Edelta(t))))+kcp)*Qs(t) + kpc*Qp(t) + kabs*Qm(t),
    Qp'(t) = -kpc*Qp(t)+kcp*Qs(t),
    Qu'(t) = (Ki+(Ks-kreabs)*CCr)*Qs(t) + k1*Qc(t),
    Qc'(t) = (-k1)*Qc(t)+V(t)*Qs(t)/(km+Qs(t)),
    SCr'(t) = -CCr/Vol*SCr(t) + k2,
    V'(t) = -alpha*V(t)*(E*gaMMa*((-k1)*Qc(t)+V(t)*Qs(t)/(km+Qs(t)))*Q50gamma/(Q50gamma+Qcgamma(t))/Qc(t)),
    Qcgamma'(t) = Qcgamma(t)*gaMMa*((-k1)*Qc(t)+V(t)*Qs(t)/(km+Qs(t)))/Qc(t),
    CS'(t) = -omega*SC(t),
    SC'(t) = omega*CS(t),
    Edelta'(t) = Edelta(t)*delta*(E*gaMMa*((-k1)*Qc(t)+V(t)*Qs(t)/(km+Qs(t)))*Q50gamma/(Q50gamma+Qcgamma(t))/Qc(t))/(Emax*Qcgamma(t)/(Q50gamma+Qcgamma(t))),
    y1(t) = SCr(t),
    y2(t) = u(t) - (Qu(t)+Qs(t)+Qp(t)+Qc(t))
)

@time println(assess_identifiability(ode))