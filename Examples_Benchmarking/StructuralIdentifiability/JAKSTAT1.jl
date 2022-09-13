#--------------------------------------------------------------------------
# Raia V et al. (2011) Dynamic mathematical modeling of IL13-induced 
# signaling in Hodgkin and primary mediastinal B-cell lymphoma allows 
# prediction of therapeutic targets. 
# Cancer Research 71(3):693{704.
#--------------------------------------------------------------------------



using Logging

using StructuralIdentifiability

logger = Logging.SimpleLogger(stdout, Logging.Info)
global_logger(logger)

ode = @ODEmodel(
    x1'(t) = -t1*x1(t)*2*u(t) - t5*x1(t) + t6*x2(t),
    x2'(t) = t5*x1(t)-t6*x2(t),
    x3'(t) = t1*2*u(t)*x1(t)-t2*x3(t)*(-x6(t)+3),
    x4'(t) = t2*x3(t)*(-x6(t)+3) - t3*x4(t),
    x5'(t) = t3*x4(t)-t4*x5(t),
    x6'(t) = -t7*x3(t)*x6(t)/(1+t13*x1(t))-t7*x4(t)*x6(t)/(1+t13*x10(t))+t8*(-x6(t)+3)*92,
    x7'(t) = -t9*x7(t)*(-x6(t)+3)+t10*(-x7(t)+165)*92,
    x8'(t) = t11*(-x7(t)+165),
    x9'(t) = -t12*2*u(t)*x9(t),
    x10'(t) = x8(t)*t14/(t15+x8(t))-t16*x10(t),
    y1(t) = x1(t)+x3(t)+x4(t),
    y2(t) = t18*(x3(t) + x4(t) + x5(t) +(1/3 -x9(t))),
    y3(t) = t19*(x4(t) + x5(t)),
    y4(t) = t20*(-x6(t)+3),
    y5(t) = t21*x8(t),
    y6(t) = t22*x8(t)*t17/t11,
    y7(t) = x10(t),
    y8(t) = -x7(t)+165
    )


@time println(assess_identifiability(ode))