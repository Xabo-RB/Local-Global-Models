using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)


ode = @ODEmodel(
    x1'(t) = -b*x1(t)*x4(t) - d*x1(t) + s,
    x2'(t) = b*q1*x1(t)*x4(t)-k1*x2(t)-w1*x2(t),
    x3'(t) = b*q2*x1(t)*x4(t)+k1*x2(t)-w2*x3(t),
    x4'(t) = -c*x4(t)+k2*x3(t),
    y1(t) = x1(t),
    y2(t) = x4(t)
)


println(find_ioequations(ode))
