#Hethcote, H., Zhien, M., & Shengbing, L. (2002).
#Effects of quarantine in six endemic models for infectious diseases. Mathematical biosciences, 180(1-2), 141-160.


using SIAN, Logging
ode = @ODEmodel(
    A'(t) = 0,
    S'(t) = A(t)-beta*S(t)*I(t)-d*S(t),
    I'(t) = (beta*S(t)-(gamma+d+delta+alpha_1))*I(t),
    Q'(t) = delta*I(t)-(epsilon+d+alpha_2)*Q(t),
    R'(t) = gamma*I(t)+epsilon*Q(t)-d*R(t),
    y1(t) = A(t),
    y2(t) = Q(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
