#P. J. Moate, R. C. Boston, T. C. Jenkins, and I. J. Lean, ‘‘Kinetics of ruminal
#lipolysis of triacylglycerol and biohydrogenation of long-chain fatty acids:
#New insights from old data,’’ J. Dairy Sci., vol. 91, no. 2, pp. 731–742,Feb. 2008

using SIAN, Logging

ode = @ODEmodel(
  x4'(t) = -k5 * x4(t) // (k6 + x4(t)),
  x5'(t) = k5 * x4(t) // (k6 + x4(t)) - k7 * x5(t) / (k8 + x5(t) + x6(t)),
  x6'(t) = k7 * x5(t) // (k8 + x5(t) + x6(t)) - k9 * x6(t) * (k10 - x6(t)) // k10,
  x7'(t) = k9 * x6(t) * (k10 - x6(t)) // k10,
  y1(t) = x4(t),
  y2(t) = x5(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p=0.99, p_mod=2^29 - 3))

using StructuralIdentifiability

ode = @ODEmodel(
  x4'(t) = -k5 * x4(t) // (k6 + x4(t)),
  x5'(t) = k5 * x4(t) // (k6 + x4(t)) - k7 * x5(t) / (k8 + x5(t) + x6(t)),
  x6'(t) = k7 * x5(t) // (k8 + x5(t) + x6(t)) - k9 * x6(t) * (k10 - x6(t)) // k10,
  x7'(t) = k9 * x6(t) * (k10 - x6(t)) // k10,
  y1(t) = x4(t),
  y2(t) = x5(t)
)
@time println(assess_identifiability(ode))