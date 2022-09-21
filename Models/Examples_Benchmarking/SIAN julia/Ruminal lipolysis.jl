#PJ Moate, RC Boston, TC Jenkins, and IJ Lean. Kinetics of ruminal lipolysis of triacylglycerol and
#biohydrogenation of long-chain fatty acids: new insights from old data. Journal of Dairy Science,
#91(2):731–742, 2008

using Logging,SIAN



ode = @ODEmodel(
    x1'(t) = -x1(t) * x5(t) / (k2 + x1(t)),
    x2'(t) = 2 * x1(t) * x5(t) / ((k2 + x1(t)) * 3) - k4 * x2(t),
    x3'(t) = k4*(x2(t))/2 - k4*x3(t),
    x4'(t) = x1(t) * x5(t) / (3 * (k2 + x1(t))) + k4 * (x2(t))/2 + k4 * x3(t),
    x5'(t) = -k3 * x5(t),
    y1(t) = x1(t),
    y2(t) = x2(t) + x3(t),
    y3(t) = x4(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))









