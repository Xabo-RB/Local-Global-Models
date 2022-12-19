#León-Triana, O., Sabir, S., Calvo, G. F., Belmonte-Beitia, J., Chulián, S., Martínez-Rubio, Á., ... & Pérez-García, V. M. (2021). CAR T cell therapy in B-cell acute lymphoblastic leukaemia: Insights from mathematical models. 
#Communications in Nonlinear Science and Numerical Simulation, 94, 105570.


using SIAN, Logging

ode = @ODEmodel(
    C'(t) = rhoc*(L(t)+B(t))*C(t) + rhob*I(t)*C(t) - C(t)/taoc, #number of cells of CAR T cells
    L'(t) = rhol*L(t) - alpha*L(t)*C(t), # number of cells of leukaemic cells
    B'(t) = I(t)/taoi - alpha*B(t)*C(t) - B(t)/taob, # number of cells of mature healthy B cells
    P'(t) = rhop*(2*ap* (1/(1 + ks*(P(t)+I(t)))) -1)*P(t) - P(t)/taop, # number of cells of CD19- haematopoietic stem cells (HSCs)
    I'(t) = rhoi*(2*ai* (1/(1 + ks*(P(t)+I(t)))) -1 )*I(t) - I(t)/taoi + P(t)/taop - alpha*beta*I(t)*C(t), # number of cells of CD19+ B cell progenitors
    y1(t) = C(t),
    y2(t) = B(t),
    y3(t) = P(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))



