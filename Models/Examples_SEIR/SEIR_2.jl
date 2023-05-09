# Assessing parameter identifiability in 
# compartmental dynamic models using a computational approach: 
# application to infectious disease transmission models 


using SIAN, Logging
ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*I(t)*S(t)/N(t),
    E'(t) =  beta*I(t)*S(t)/N(t)-k*E(t),
    I'(t) = k*E(t)-gamma*I(t),
    R'(t) = gamma*I(t),
    C'(t) = k*E(t),
    y1(t) = kk*C(t),
    y2(t) = N(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

