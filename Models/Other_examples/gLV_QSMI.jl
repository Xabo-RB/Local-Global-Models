#Brunner, J. D., & Chia, N. (2019). Metabolite-mediated modelling of microbial community dynamics captures emergent behaviour more effectively than species–species modelling.
# Journal of the Royal Society Interface, 16(159), 20190423.


using SIAN, Logging
ode = @ODEmodel(
    x1'(t) = k11*x1(t)*w1(t) - d1*x1(t) + phi12*x1(t)*w2(t),
    x2'(t) = k21*x2(t)*w1(t) - d2*x(2),
    w1'(t) = f1 - dd1*w1(t) - k11*x1(t)*w1(t) - k21*x2(t)*w1(t),
    w2'(t) = k21*x2(t)*w1(t) - dd2*w2(t) - k12*x1(t)*w2(t),
    #yy1'(t) = k11*x1(t)*y1(t) - ddd1*yy1(t) - kk21*x2(t)*yy1(t),
    #yy2'(t) = kk21*x2(t)*yy1(t) - ddd2*yy2(t) - kk32*x3(t)*yy2(t),
    y1(t) = x1(t),
    y2(t) = x2(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
