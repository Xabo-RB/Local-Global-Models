#Garcia, V., Bonhoeffer, S., & Fu, F. (2020). Cancer-induced immunosuppression can enable effectiveness of immunotherapy through
#bistability generation: A mathematical and computational examination. Journal of theoretical biology, 492, 110185.


using SIAN, Logging

ode = @ODEmodel(
    x'(t) = (r1 + b1*y(t))*x(t)*(1-b1*x(t))- d1*x(t)*z(t)/(m1+w(t)), #pancreatic cancer cell population
    y'(t) =  (r2 + b2*w(t)/(k2+w(t)))*y(t)*(1-b2*y(t))-mu2*y(t), #pancreatic stellate cell population
    y1(t) = z(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

