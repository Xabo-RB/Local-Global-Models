read "C:/Users/Even/OneDrive - Universidade de Vigo/Escritorio/Softwares Benchmarking/SIAN-master/IdentifiabilityODE.mpl";


sigma := [
  diff(S(t), t) = -betta*S(t)*II(t),
  diff(E(t), t) =  betta*S(t)*II(t)-v*E(t),
  diff(II(t), t) = v*E(t)-pssi*II(t)-(1-pssi)*gammma*II(t),
  diff(R(t), t) = gammma*Q(t)+(1-pssi)*gammma*II(t),
  diff(Q(t), t) = -gammma*Q(t)+pssi*II(t),
  y1(t) = Q(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), infolevel = 2):