#Lan K Nguyen, Andrea Degasperi, Philip Cotter, and Boris N Kholodenko. Dyvipac: an integrated
#analysis and visualisation framework to probe multi-dimensional biological networks. Scientific reports,
#5(1):1–17, 2015


using Logging,SIAN



ode = @ODEmodel(
    x1'(t) = k1*(s1t-x1(t))*(KK1^n1)/(KK1^n1+x3(t)^n1)-k2*x1(t),
    x2'(t) = k3*(s2t-x2(t))*x1(t)*(1+(alpha*x3(t)^n2)/(KK2^n2+x3(t)^n2))-k4*x2(t),
    x3'(t) = k5*(s3t-x3(t))*x2(t)-k6*x3(t),
    y1(t) = x1(t),
    y2(t) = x2(t),
    y3(t) = x3(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1))












