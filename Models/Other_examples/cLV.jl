#Joseph, T. A., Shenhav, L., Xavier, J. B., Halperin, E., & Pe’er, I. (2020). Compositional Lotka-Volterra describes microbial dynamics in the simplex.
# PLoS computational biology, 16(5), e1007917.


using SIAN, Logging
ode = @ODEmodel(
    pi1'(t) = pi1(t)*( (g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    pi2'(t) = pi2(t)*( (g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    pi3'(t) = pi3(t)*( (g3+A31*pi1(t)+A32*pi2(t)+A33*pi3(t)+B31*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    y1(t) = pi1(t),
    y2(t) = pi2(t),
    y3(t) = pi3(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#Model without external inputs
#2
using SIAN, Logging
ode = @ODEmodel(
    pi1'(t) = pi1(t)*( (g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t))),
    pi2'(t) = pi2(t)*( (g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t))),
    pi3'(t) = pi3(t)*( (g3+A31*pi1(t)+A32*pi2(t)+A33*pi3(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t))),
    y1(t) = pi1(t),
    y2(t) = pi2(t),
    y3(t) = pi3(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#Model without external inputs 1 less output
#3
using SIAN, Logging
ode = @ODEmodel(
    pi1'(t) = pi1(t)*( (g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t))),
    pi2'(t) = pi2(t)*( (g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t))),
    pi3'(t) = pi3(t)*( (g3+A31*pi1(t)+A32*pi2(t)+A33*pi3(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t))),
    y1(t) = pi1(t),
    y2(t) = pi2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


# 1 less output
#4
using SIAN, Logging
ode = @ODEmodel(
    pi1'(t) = pi1(t)*( (g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    pi2'(t) = pi2(t)*( (g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    pi3'(t) = pi3(t)*( (g3+A31*pi1(t)+A32*pi2(t)+A33*pi3(t)+B31*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    y1(t) = pi1(t),
    y2(t) = pi2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

# 1 less output
using SIAN, Logging
ode = @ODEmodel(
    pi1'(t) = pi1(t)*( (g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    pi2'(t) = pi2(t)*( (g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    pi3'(t) = pi3(t)*( (g3+A31*pi1(t)+A32*pi2(t)+A33*pi3(t)+B31*u1(t)) - pi1(t)*(g1+A11*pi1(t)+A12*pi2(t)+A13*pi3(t)+B11*u1(t)) + pi2(t)*(g2+A21*pi1(t)+A22*pi2(t)+A23*pi3(t)+B21*u1(t))),
    y1(t) = pi1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
