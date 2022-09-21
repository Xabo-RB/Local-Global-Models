#Alan S Perelson and Patrick W Nelson. Mathematical analysis
#of hiv-1 dynamics in vivo. SIAM review,
#41(1):3–44, 1999

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)


ode = @ODEmodel(
    x1'(t) = lambda-rho*x1(t)-eta(t)*x1(t)*x3(t),
    x2'(t) = eta(t)*x1(t)*x3(t)-delta*x2(t),
    x3'(t) = N*delta*x2(t)-c*x3(t),
    y1(t) = x3(t),
    y2(t) = x1(t) + x2(t),
)

@time println(assess_identifiability(ode))