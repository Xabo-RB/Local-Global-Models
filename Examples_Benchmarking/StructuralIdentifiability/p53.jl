#Joseph Distefano. Dynamic systems biology modeling and simulation. Academic Press, 2015

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = (p1*x4(t)) - (p3*x1(t)) - p4*( (x1(t)^2/(p5+x1(t)))* (1+ (p6*u1(t)/(p7+u1(t))) ) ),
    x2'(t) = p8 - (p9*x2(t))- p10*( (x1(t)*x2(t)/(p11+x2(t))) * (1+(p12*u1(t)/(p13+u1(t))) ) ),
    x3'(t) = p14 - (p15*x3(t)) - p16*x1(t)*x3(t)*(1-p18*u1(t))/(p17+x3(t)),
    x4'(t) = p20 - p21*(1-p24)*(1-p25)/((p22^4)+1) - (p20*x4(t)) + (p21*(x3(t)^4)) + (1+p23*u1(t))*(1-p24*x1(t))*(1-p25*x2(t))/(p22^4+x3(t)^4),
    y1(t) = x1(t),
    y2(t) = x2(t),
    y3(t) = x3(t),
    y4(t) = x4(t)
)

@time println(assess_identifiability(ode))












