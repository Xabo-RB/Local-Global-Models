using SIAN, Logging
ode = @ODEmodel(
    Cr'(t) = q*(u1(t)-Cr(t))-0.072*(exp(gamma*Tr(t)/(gamma+Tr(t)/20)))*Cr(t)+d,
    Tr'(t) =  q*(u2(t)-Tr(t))+0.576*(exp(gamma*Tr(t)/(gamma+Tr(t)/20)))*Cr(t)-0.3*(Tr(t)-Tc(t)),
    Tc'(t) = delta1*qc*(u3(t)-Tc(t))+3.0*(Tr(t)-Tc(t)),
    
    y1(t) = Tc(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

##NO FUNCIONA TIENE EXPONENCIALES