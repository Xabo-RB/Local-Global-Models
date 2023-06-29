using SIAN, Logging

ode = @ODEmodel(
  mRNA'(t) = -d*mRNA,
  GFP'(t) = kTL*mRNA*u1(t) - b*GFP*u2(t),
  y1(t) = GFP(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


using StructuralIdentifiability

ode = @ODEmodel(
    mRNA'(t) = -d*mRNA,
    GFP'(t) = kTL*mRNA*u1(t) - b*GFP*u2(t),
    y1(t) = GFP(t)
)

@time println(assess_identifiability(ode))
