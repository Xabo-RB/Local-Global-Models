# Example 6.3 from the paper "Global Identifiability of Differential Models", taken from
# Balsa-Canto, E., Alonso, A. A., Banga, J. R., 
# An iterative identification procedure for dynamic modeling of biochemical networks
st := time()
read "C:/Users/Even/OneDrive - Universidade de Vigo/Escritorio/Softwares Benchmarking/SIAN-master/IdentifiabilityODE.mpl";


sigma := [
  diff(x1(t), t) = lambda-rho*x1(t)-eta(t)*x1(t)*x3(t),
  diff(x2(t), t) = eta(t)*x1(t)*x3(t)-delta*x2(t),
  diff(x3(t), t) = N*delta*x2(t)-c*x3(t),
  y1(t) = x3(t),
  y2(t) = x1(t) + x2(t),
];

IdentifiabilityODE(sigma, GetParameters(sigma));
time() - st;