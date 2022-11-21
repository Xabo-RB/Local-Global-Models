#Sousa, J., & Carneiro, J. (2000). A mathematical analysis of TCR serial triggering and down‐regulation. 
#European Journal of Immunology, 30(11), 3219-3227.


using SIAN, Logging

ode = @ODEmodel(
    T'(t) = - k * (T(t)^2)*(L^2),
    y1(t) = T(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
