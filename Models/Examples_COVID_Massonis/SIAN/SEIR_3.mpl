read "C:/Users/Even/OneDrive - Universidade de Vigo/Escritorio/Softwares Benchmarking/SIAN-master/IdentifiabilityODE.mpl";


sigma := [
  diff(S(t), t) = betta*II(t)*S(t)/N(t)-allpha*S(t),
  diff(E(t), t) = betta*II(t)*S(t)/N(t)-gammma*E(t),
  diff(II(t), t) = gammma*E(t)-dellta*II(t),
  diff(R(t), t) = dellta*II(t)-allpha*Q(t)-Kk*Q(t),
  diff(Q(t), t) = allpha*Q(t),
  diff(Dd(t), t) = Kk*Q(t),
  diff(P(t), t) = allpha*S(t),
  diff(N(t), t) = 0,
  y1(t) = Q(t),
  y2(t) = N(t),
  y3(t) = Dd(t),
  y4(t) = R(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), infolevel = 2):