#Hu, X., Ke, G., & Jang, S. R. J. (2019). 
#Modeling pancreatic cancer dynamics with immunotherapy. Bulletin of mathematical biology, 81(6), 1885-1915.

#NO FUNCIONA
using SIAN, Logging

ode = @ODEmodel(
    x'(t) = (r1 + b1*y(t))*x(t)*(1-b1*x(t))- d1*x(t)*z(t)/(m1+w(t)), #pancreatic cancer cell population
    y'(t) =  (r2 + b2*w(t)/(k2+w(t)))*y(t)*(1-b2*y(t))-mu2*y(t), #pancreatic stellate cell population
    z'(t) = b3*z(t)*v(t)/((k3+v(t))*(m3+w(t))) - mu3*z(t) +r3, #effector cells, including CD8+T cells and NK cells
    w'(t) = b4*x(t)*z(t)/(k4+x(t)) - mu4*w(t) + r4*x(t)*y(t)/(m4+v(t)), #concentration of tumor promoting cytokines, including TGF-beta and IL-6
    v'(t) = b5*x(t)*z(t)/(k5 + x(t))- mu5*v(t), #concentration of tumor suppressing cytokines, including INF-gamma and IL-2
    y1(t) = z(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#NO FUNCIONA
using SIAN, Logging

ode = @ODEmodel(
    x'(t) = (r1 + b1*y(t))*x(t)*(1-b1*x(t))- d1*x(t)*z(t)/(m1+w(t)), #pancreatic cancer cell population
    y'(t) =  (r2 + b2*w(t)/(k2+w(t)))*y(t)*(1-b2*y(t))-mu2*y(t), #pancreatic stellate cell population
    z'(t) = b3*z(t)*v(t)/((k3+v(t))*(m3+w(t))) - mu3*z(t) +r3, #effector cells, including CD8+T cells and NK cells
    w'(t) = b4*x(t)*z(t)/(k4+x(t)) - mu4*w(t) + r4*x(t)*y(t)/(m4+v(t)), #concentration of tumor promoting cytokines, including TGF-beta and IL-6
    v'(t) = b5*x(t)*z(t)/(k5 + x(t))- mu5*v(t), #concentration of tumor suppressing cytokines, including INF-gamma and IL-2
    y1(t) = z(t),
    y2(t) = y(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#NO FUNCIONA
using SIAN, Logging

ode = @ODEmodel(
    x'(t) = (r1 + b1*y(t))*x(t)*(1-b1*x(t))- d1*x(t)*z(t)/(m1+w(t)), #pancreatic cancer cell population
    y'(t) =  (r2 + b2*w(t)/(k2+w(t)))*y(t)*(1-b2*y(t))-mu2*y(t), #pancreatic stellate cell population
    z'(t) = b3*z(t)*v(t)/((k3+v(t))*(m3+w(t))) - mu3*z(t) +r3, #effector cells, including CD8+T cells and NK cells
    w'(t) = b4*x(t)*z(t)/(k4+x(t)) - mu4*w(t) + r4*x(t)*y(t)/(m4+v(t)), #concentration of tumor promoting cytokines, including TGF-beta and IL-6
    v'(t) = b5*x(t)*z(t)/(k5 + x(t))- mu5*v(t), #concentration of tumor suppressing cytokines, including INF-gamma and IL-2
    y1(t) = z(t),
    y2(t) = y(t) + x(t),
    y3(t) = v(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#OUTPUT 2
using SIAN, Logging

ode = @ODEmodel(
    x'(t) = (r1 + b1*y(t))*x(t)*(1-b1*x(t))- d1*x(t)*z(t)/(m1+w(t)), #pancreatic cancer cell population
    y'(t) =  (r2 + b2*w(t)/(k2+w(t)))*y(t)*(1-b2*y(t))-mu2*y(t), #pancreatic stellate cell population
    z'(t) = b3*z(t)*v(t)/((k3+v(t))*(m3+w(t))) - mu3*z(t) +r3, #effector cells, including CD8+T cells and NK cells
    w'(t) = b4*x(t)*z(t)/(k4+x(t)) - mu4*w(t) + r4*x(t)*y(t)/(m4+v(t)), #concentration of tumor promoting cytokines, including TGF-beta and IL-6
    v'(t) = b5*x(t)*z(t)/(k5 + x(t))- mu5*v(t), #concentration of tumor suppressing cytokines, including INF-gamma and IL-2
    y1(t) = z(t),
    y2(t) = y(t),
    y3(t) = x(t),
    y4(t) = v(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10))

