#Alejandro F Villaverde, Neil D Evans, Michael J Chappell, and Julio R Banga. Input-dependent
#structural identifiability of nonlinear systems. IEEE Control Systems Letters, 3(2):272–277, 2018

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = -(k12 + k1e)*x1(t) + k21*x2(t) + b*u(t),
    x2'(t) = k12 * x1(t) - k21 * x2(t),
    y1(t) = x1
)

@time println(assess_identifiability(ode))

#u constant
using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = -(k12 + k1e)*x1(t) + k21*x2(t) + b*u,
    x2'(t) = k12 * x1(t) - k21 * x2(t),
    y1(t) = x1
)

@time println(assess_identifiability(ode))
