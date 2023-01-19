st := time()
read "C:/Users/Even/OneDrive - Universidade de Vigo/Escritorio/Softwares Benchmarking/SIAN-master/IdentifiabilityODE.mpl";


sigma := [
  diff(S(t), t) = -beta*S(t)*I(t),
  diff(E(t), t) =  beta*S(t)*I(t)-v*E(t),
  diff(I(t), t) = v*E(t)-psi*I(t)-(1-psi)*gamma*I(t),
  diff(R(t), t) = gamma*Q(t)+(1-psi)*gamma*I(t),
  diff(Q(t), t) = -gamma*Q(t)+psi*I(t),
  y1(t) = Q(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), infolevel = 2);
time() - st;