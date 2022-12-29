#Qiao, Y. H., Xu, C. Q., Zeng, Y. J., Xu, X. H., Zhao, H., & Xu, H. (2004). The kinetic model and simulation of blood coagulation—the kinetic influence of activated protein C. 
#Medical engineering & physics, 26(4), 341-347.

using SIAN, Logging
ode = @ODEmodel(
    IXa'(t) = k1*beta - h1*IXa(t),
    VIIIa'(t) = k2*IIa(t) + k3*APC(t)*VIIIa(t)/(b1+VIIIa(t)) - h2*VIIIa(t),
    Xa'(t) = k5*IXa(t)*VIIIa(t)/(b2+VIIIa(t)) - h3*Xa(t),
    Va'(t) = k6*IIa(t) - k7*APC(t)*Va(t)/(b3+Va(t)) - h4*Va(t),
    APC'(t) = k8*IIa(t) - h5*APC(t),
    IIa'(t) = k9*Xa(t)*Va(t)/(b4+Va(t))-h6*IIa(t),
    y1(t) = Xa(t),
    y2(t) = IIa(t),
    y3(t) = IXa(t)
)
@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
