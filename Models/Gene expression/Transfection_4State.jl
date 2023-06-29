using SIAN, Logging

ode = @ODEmodel(
  mRNA'(t) = -d1*mRNA(t)-d2*mRNA(t)*enz(t),
  GFP'(t) = kTL*mRNA(t)-b*GFP(t),
  enz'(t) = d3*mRNAenz(t)-d2*mRNA(t)*enz(t),
  mRNAenz'(t) = -d3*mRNAenz(t)+d2*mRNA(t)*enz(t),
  y1(t) = GFP(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
