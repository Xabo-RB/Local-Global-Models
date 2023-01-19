read "C:/Users/Even/OneDrive - Universidade de Vigo/Escritorio/Softwares Benchmarking/SIAN-master/IdentifiabilityODE.mpl";


sigma := [
  diff(S(t), t) = -betta*II(t)*S(t)/N(t),
  diff(E(t), t) =  betta*II(t)*S(t)/N(t)-k*E(t),
  diff(II(t), t) = k*E(t)-gammma*II(t),
  diff(R(t), t) = gammma*II(t),
  diff(C(t), t) = k*E(t),
  diff(N(t), t) = 0,
  y1(t) = kk*C(t),
  y2(t) = N(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), infolevel = 2):