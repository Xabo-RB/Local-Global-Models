# Li, G., Leung, C. Y., Wardi, Y., Debarbieux, L., & Weitz, J. S. (2020). Optimizing the timing and composition of therapeutic phage cocktails: 
#a control-theoretic approach. Bulletin of mathematical biology, 82(6), 1-29.

using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = r1*x1(t)*(1-((x1(t)+x2(t))/kCD))*(1-m)-kPD*x1(t)*Y*x3(t)/(1+x3(t))-((E1*x4(t)*x1(t))/(1+x1(t)+x2(t))),
    x2'(t) = r2*x2(t)*(1-((x1(t)+x2(t))/kCD))+r1*x1(t)*(1-((x1(t)+x2(t))/kCD))*m-kPD*x2(t)*Y*x5(t)/(1+x5(t))-((E1*x4(t)*x2(t))/(1+x1(t)+x2(t))),
    x3'(t) = beta*x1(t)*Y*x3(t)/(1+x3(t))-Y*x1(t)*x3(t)-v*x3(t)+q*u1(t),
    x4'(t) = alpha*x4(t)*(1-x4(t))*((x1(t)+x2(t))/(x1(t)+x2(t)+kND)),
    x5'(t) = beta*x2(t)*Y*x5(t)/(1+x5(t))-Y*x2(t)*x5(t)-v*x5(t)+q*u2(t),
    y1(t) = x1(t),
    y2(t) = x2(t),
    y3(t) = x3(t),
    y4(t) = x4(t),
    y5(t) = x5(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
