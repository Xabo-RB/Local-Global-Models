#Garcia, V., Bonhoeffer, S., & Fu, F. (2020). Cancer-induced immunosuppression can enable effectiveness of immunotherapy through
#bistability generation: A mathematical and computational examination. Journal of theoretical biology, 492, 110185.


using SIAN, Logging

ode = @ODEmodel(
    T'(t) = a*T(t)*(1-b*T(t))-k*T(t)*E(t), #Tumor cells
    E'(t) =  sigma - d*E(t)+m*E(t)*T(t), # effector cells
    y1(t) = E(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10))

