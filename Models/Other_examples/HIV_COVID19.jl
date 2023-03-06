#MODEL SELECTION AND IDENTIFIABILITY ANALYSIS OF HIVE AND SARS COV2 CO INFECTION MODEL WITH DRUG THERAPY

using SIAN, Logging
ode = @ODEmodel(
    T'(t) = lambdaH - betaH*T(t)*Vh(t)*hrt + delta*T(t)*C(t) - dT*T(t),
    I'(t) = betaH*T(t)*Vh(t)*hrt - dI*I(t),
    Vh'(t) = sigmaH*I(t)-cH*Vh(t),
    P'(t) = lambdaS - betaM*P(t)*(Pm(t)+Vm(t)) - dP*P(t),
    Pm'(t) = betaM*P(t)*(Pm(t)+Vm(t))-dS*Pm(t),
    Vm'(t) = sigmaM*Pm(t)-cm*Vm(t),
    C'(t) = mu*T(t)*(Pm(t)+Vm(t)) - alpha*C(t) - dC*C(t),
    y1(t) = T(t)+I(t),
    y2(t) = Vh(t),
    y3(t) = Vm(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#[ Info: === Summary ===
#[ Info: Globally identifiable parameters:                 [T(0), I(0), Vh(0), Vm(0), lambdaH, dT, dI, sigmaH, cH, dP]
#[ Info: Locally but not globally identifiable parameters: [Pm(0), betaM, lambdaS, dS, sigmaM, cm, P(0)]

#[ Info: Not identifiable parameters:                      [C(0), betaH, hrt, delta, mu, alpha, dC]

#Los seis parámetros globalmente identificables se corresponde con aquellos que en el paper considera que son identificable. Los
#SLI no aparecen.

#ORDENADOR GRANDE
using StructuralIdentifiability

ode = @ODEmodel(
    T'(t) = lambdaH - betaH*T(t)*Vh(t)*hrt + delta*T(t)*C(t) - dT*T(t),
    I'(t) = betaH*T(t)*Vh(t)*hrt - dI*I(t),
    Vh'(t) = sigmaH*I(t)-cH*Vh(t),
    P'(t) = lambdaS - betaM*P(t)*(Pm(t)+Vm(t)) - dP*P(t),
    Pm'(t) = betaM*P(t)*(Pm(t)+Vm(t))-dS*Pm(t),
    Vm'(t) = sigmaM*Pm(t)-cm*Vm(t),
    C'(t) = mu*T(t)*(Pm(t)+Vm(t)) - alpha*C(t) - dC*C(t),
    y1(t) = Vh(t)
)

println(find_ioequations(ode))


#1 sólo output / ORDENADOR GRANDE
using SIAN, Logging
ode = @ODEmodel(
    T'(t) = lambdaH - betaH*T(t)*Vh(t)*hrt + delta*T(t)*C(t) - dT*T(t),
    I'(t) = betaH*T(t)*Vh(t)*hrt - dI*I(t),
    Vh'(t) = sigmaH*I(t)-cH*Vh(t),
    P'(t) = lambdaS - betaM*P(t)*(Pm(t)+Vm(t)) - dP*P(t),
    Pm'(t) = betaM*P(t)*(Pm(t)+Vm(t))-dS*Pm(t),
    Vm'(t) = sigmaM*Pm(t)-cm*Vm(t),
    C'(t) = mu*T(t)*(Pm(t)+Vm(t)) - alpha*C(t) - dC*C(t),
    y1(t) = Vh(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#Hrt(t), mismo resultado como un estado no medible y de valor constante, no tiene derivadas
using SIAN, Logging
ode = @ODEmodel(
    T'(t) = lambdaH - betaH*T(t)*Vh(t)*hrt(t) + delta*T(t)*C(t) - dT*T(t),
    I'(t) = betaH*T(t)*Vh(t)*hrt(t) - dI*I(t),
    Vh'(t) = sigmaH*I(t)-cH*Vh(t),
    P'(t) = lambdaS - betaM*P(t)*(Pm(t)+Vm(t)) - dP*P(t),
    Pm'(t) = betaM*P(t)*(Pm(t)+Vm(t))-dS*Pm(t),
    Vm'(t) = sigmaM*Pm(t)-cm*Vm(t),
    C'(t) = mu*T(t)*(Pm(t)+Vm(t)) - alpha*C(t) - dC*C(t),
    hrt'(t) = 0,
    y1(t) = T(t)+I(t),
    y2(t) = Vh(t),
    y3(t) = Vm(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#1 derivada -> betaH y hrt pasan a ser SGI cuando Hrt(t) y no es de valor constante. Una entrada desconocida/ un estado no medible?
using SIAN, Logging
ode = @ODEmodel(
    T'(t) = lambdaH - betaH*T(t)*Vh(t)*hrt(t) + delta*T(t)*C(t) - dT*T(t),
    I'(t) = betaH*T(t)*Vh(t)*hrt(t) - dI*I(t),
    Vh'(t) = sigmaH*I(t)-cH*Vh(t),
    P'(t) = lambdaS - betaM*P(t)*(Pm(t)+Vm(t)) - dP*P(t),
    Pm'(t) = betaM*P(t)*(Pm(t)+Vm(t))-dS*Pm(t),
    Vm'(t) = sigmaM*Pm(t)-cm*Vm(t),
    C'(t) = mu*T(t)*(Pm(t)+Vm(t)) - alpha*C(t) - dC*C(t),
    hrt'(t) = 1,
    y1(t) = T(t)+I(t),
    y2(t) = Vh(t),
    y3(t) = Vm(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#[ Info: === Summary ===
#[ Info: Globally identifiable parameters:                 [lambdaH, dI, cH, hrt, I, Vm, betaH, dT, sigmaH, dP, T, Vh]
#[ Info: Locally but not globally identifiable parameters: [Pm, betaM, lambdaS, dS, sigmaM, cm, P]
#[ Info: Not identifiable parameters:                      [delta, mu, alpha, dC, C]