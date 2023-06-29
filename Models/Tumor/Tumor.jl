#Maria Pia Saccomani, Stefania Audoly, Giuseppina Bellu, and Leontina D’Angiò. Examples of testing
#global identifiability of biological and biomedical models with the daisy software. Computers in Biology
#and Medicine, 40(4):402–407, 2010

using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = -(k3+k7)*x1(t) +k4*x2(t),
    x2'(t) = k3*x1(t)-(k4+a*k5+b*d*k5)*x2(t)+k6*x3(t)+k6*x4(t)+k5*x2(t)*x3(t)+k5*x2(t)*x4(t),
    x3'(t) = a*k5*x2(t)-k6*x3(t)-k5*x2(t)*x4(t),
    x4'(t) = b*d*k5*x2(t)-k6*x4(t)-k5*x2(t)*x4(t),
    x5'(t) = k7*x1(t),
    y1(t) = x5(t)
)

@time println(assess_identifiability(ode))





