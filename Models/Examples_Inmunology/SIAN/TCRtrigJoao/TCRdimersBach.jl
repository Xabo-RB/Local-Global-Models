#Bachmann, M. F., Salzmann, M., Oxenius, A., & Ohashi, P. S. (1998). Formation of TCR dimers/trimers as a crucial step for T cell activation.
#European journal of immunology, 28(8), 2571-2579.


using SIAN, Logging

ode = @ODEmodel(
    T'(t) = - keff * (T(t)^2)*((kon/koff)^2),
    y1(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
