using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = -e*x1(t) - a*x1(t)*(p-x2(t)) + d*x2(t),
    x2'(t) = a*x1(t)*(p-x2(t)) - d*x2(t) ,
    y1(t) = c*(x1(t)+x2(t))
)

@time println(assess_identifiability(ode))