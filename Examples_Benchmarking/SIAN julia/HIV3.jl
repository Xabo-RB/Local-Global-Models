# Wodarz, D., Nowak, M.
# Specific therapy regimes could lead to long-term immunological control of HIV
# https://doi.org/10.1073/pnas.96.25.14464
# Page 1

using Logging, SIAN


ode = @ODEmodel(
    x'(t) = lm - d * x(t) - beta * x(t) * v(t),
    y'(t) = beta * x(t) * v(t) - a * y(t),
    v'(t) = k * y(t) - u * v(t),
    w'(t) = c * x(t) * y(t) * w(t) - c * q * y(t) * w(t) - b * w(t),
    z'(t) = c * q * y(t) * w(t) - h * z(t),
    y1(t) = w(t),
    y2(t) = z(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, infolevel = 10, nthrds = 1))