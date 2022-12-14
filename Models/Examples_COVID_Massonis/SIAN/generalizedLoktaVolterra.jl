# Remien, C. H., Eckwright, M. J., & Ridenhour, B. J. (2021). Structural identifiability of the 
#generalized Lotka–Volterra model for microbiome studies. Royal Society Open Science, 8(7), 201378.

using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = r1*x1(t)+beta11*x1(t)^2+beta12*x1(t)*x2(t),
    x2'(t) = r2*x2(t)+beta21*x1(t)*x2(t)+beta22*x2(t)^2,
    y1(t) = x1(t),
    y2(t) = x2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
