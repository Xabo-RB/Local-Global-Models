#León-Triana, O., Sabir, S., Calvo, G. F., Belmonte-Beitia, J., Chulián, S., Martínez-Rubio, Á., ... & Pérez-García, V. M. (2021). CAR T cell therapy in B-cell acute lymphoblastic leukaemia: Insights from mathematical models. 
#Communications in Nonlinear Science and Numerical Simulation, 94, 105570.


using SIAN, Logging

ode = @ODEmodel(
    x'(t) = (r1 + b1*y(t))*x(t)*(1-b1*x(t))- d1*x(t)*z(t)/(m1+w(t)), #
    y'(t) =  (r2 + b2*w(t)/(k2+w(t)))*y(t)*(1-b2*y(t))-mu2*y(t), #
    z'(t) = b3*z(t)*v(t)/((k3+v(t))*(m3+w(t))) - mu3*z(t) +r3, #
    w'(t) = b4*x(t)*z(t)/(k4+x(t)) - mu4*w(t) + r4*x(t)*y(t)/(m4+v(t)), #
    v'(t) = b5*x(t)*z(t)/(k5 + x(t))- mu5*v(t), #
    y1(t) = z(t),
    y2(t) = y(t),
    y3(t) = x(t),
    y4(t) = v(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))



