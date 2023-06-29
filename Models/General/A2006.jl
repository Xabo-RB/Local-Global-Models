using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x10'(t) = x21(t) - x11(t),
    x11'(t) = (x10(t)-x11(t))/r,
    x20'(t) = - x21(t),
    x21'(t) = (x20(t)-x21(t))/r,
    y1(t) = x10(t)
)

@time println(assess_identifiability(ode))
