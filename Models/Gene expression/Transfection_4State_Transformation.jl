using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -x1(t)*(d1 + d2_times_mRNA0*x3(t)),
  GFP'(t) = kTL_times_mRNA0*x1(t) - GFP*b,
  x3'(t) = d3*enz0_div_mRNA0 - d3*x3(t) - d2_times_mRNA0*x3(t)*x1(t),
  y1(t) = GFP(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
