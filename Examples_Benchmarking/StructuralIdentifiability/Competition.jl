

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = x1*(a1-b11*log(x1(t))-b12*log(x2(t))),
    x2'(t) = x2(t)*(a2-b21*log(x1(t))-b22*log(x2(t))),
    y1(t) = x1
)

@time println(assess_identifiability(ode))


