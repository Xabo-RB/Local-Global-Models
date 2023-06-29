#Atilla Raksanyi. Utilisation du calcul formel pour l’étude des systèmes d’équations polynomiales (ap-
#plications en modélisation). PhD thesis, Paris 9, 1986


using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)


ode = @ODEmodel(
    x1'(t) = u1(t)-(k1+k2)*x1(t),
    x2'(t) = k1*x1(t)-(k3+k6+k7)*x2(t)+k5*x4(t),
    x3'(t) = k2*x1(t)+k3*x2(t)-k4*x3(t),
    x4'(t) = k6*x2(t)-k5*x4(t),
    y1(t) = s2*x2(t),
    y2(t) = s3*x3(t)
)

@time println(assess_identifiability(ode))










