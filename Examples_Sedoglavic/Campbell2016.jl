using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = -p_3* x1(t) - p_4* x1(t)^2 + p_1 *x2(t) + p_2 *x3(t) - p_7 *x1(t) *u(t),
    x2'(t) = p_5* p_6* x1(t) - p_3* x2(t) - p_4* x1(t)* x2(t) - p_5* x2(t) - p_7* x2(t)* u(t),
    x3'(t) = 2* p_5* p_6* x2(t) - p_3* x3(t) - p_4* x1(t)* x3(t) - 2* p_5* x3(t) - p_7* x3(t) *u(t),
    y1(t) = p_7* x1(t)* u(t),
    y2(t) = p_7* x2(t)* u(t),
    y3(t) = p_7* x3(t)* u(t)
)

@time println(assess_identifiability(ode))
