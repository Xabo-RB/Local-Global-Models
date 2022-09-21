#Carsten Conradi and Anne Shiu. Dynamics of posttranslational modification systems: Recent progress
#and future directions. Biophysical journal, 114(3):507–515, 2018

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)


ode = @ODEmodel(
    x5'(t) = k5*x6(t) + k4*x6(t) - k6*x5(t)*x3(t),
    x6'(t) = -k5*x6(t) - k4*x6(t) + k6*x5(t)*x3(t),
    x4'(t) = -k3*x4(t) - k2*x4(t) + k1*x1(t)*x2(t),
    x2'(t) = k3*x4(t) + k2*x4(t) + k1*x1(t)*x2(t),
    x1'(t) = k4*x6(t) + k2*x4(t) - k1*x1(t)*x2(t),
    x3'(t) = k5*x6(t) + k3*x4(t) - k6*x5(t)*x3(t),
    y1(t) = x3(t),
    y2(t) = x2(t)
)


@time println(assess_identifiability(ode))