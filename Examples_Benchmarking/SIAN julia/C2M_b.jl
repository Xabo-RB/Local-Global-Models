#Alejandro F Villaverde, Neil D Evans, Michael J Chappell, and Julio R Banga. Input-dependent
#structural identifiability of nonlinear systems. IEEE Control Systems Letters, 3(2):272–277, 2018

using Logging, SIAN

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

#u constant

ode = @ODEmodel(
    x1'(t) = -(k12 + k1e)*x1(t) + k21*x2(t) + b*u,
    x2'(t) = k12 * x1(t) - k21 * x2(t),
    y1(t) = x1
)

@time res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)
