#de Pillis, L. G., Fister, K. R., Gu, W., Head, T., Maples, K., Neal, T., ... & Kozai, K. (2008). Optimal control of mixed immunotherapy and chemotherapy of tumors.
# Journal of Biological systems, 16(01), 51-80.

#No funciona
using SIAN, Logging

ode = @ODEmodel(
    T'(t) = a*T(t)*(1-b*T(t))-c1*N(t)*T(t)-D(t)*T(t)-KT*M(t)*T(t), #tumor cells
    L'(t) = -m*L(t) - q*L(t)*T(t)-ucte*L(t)^2+r2*C(t)*T(t) + pI*L(t)*I(t)/(gI+I(t)) + u1(t) - KL*M(t)*L(t), # tumor-specific effector cells, T-celss
    N'(t) = alpha1 - f*N(t)+ g*T(t)*N(t)/(h + T(t)) - p*N(t)*T(t) - KN*M(t)*N(t), # non-specific effector cells, NK cells
    C'(t) = alpha2 - beta*C(t) - KC*M(t)*C(t), #circulating lymphocytes
    I'(t) = pt*T(t)*L(t)/(gt + T(t)) +  w*L(t)*I(t) - muI*I(t) + u2(t), # IL-2, VI = u2 aplicación directa, terapia de IL2
    M'(t) = - gamma*M(t) + u1(t), #chemotherapy drug, terapia/aplicación de quimio, u1 = VM
    y1(t) = L(t),
    y2(t) = N(t),
    y3(t) = M(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))



#NO FUNCIONA
using SIAN, Logging

ode = @ODEmodel(
    T'(t) = a*T(t)*(1-b*T(t))-c1*N(t)*T(t)-D(t)*T(t)-KT*M(t)*T(t), #tumor cells
    L'(t) = -m*L(t) - q*L(t)*T(t)-ucte*L(t)^2+r2*C(t)*T(t) + pI*L(t)*I(t)/(gI+I(t)) + u1(t) - KL*M(t)*L(t), # tumor-specific effector cells, T-celss
    N'(t) = alpha1 - f*N(t)+ g*T(t)*N(t)/(h + T(t)) - p*N(t)*T(t) - KN*M(t)*N(t), # non-specific effector cells, NK cells
    C'(t) = alpha2 - beta*C(t) - KC*M(t)*C(t), #circulating lymphocytes
    I'(t) = pt*T(t)*L(t)/(gt + T(t)) +  w*L(t)*I(t) - muI*I(t) + u2(t), # IL-2, VI = u2 aplicación directa, terapia de IL2
    M'(t) = - gamma*M(t) + u1(t), #chemotherapy drug, terapia/aplicación de quimio, u1 = VM
    y1(t) = L(t) + N(t) + C(t),
    y2(t) = M(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#NO FUNCIONA
using SIAN, Logging

ode = @ODEmodel(
    T'(t) = a*T(t)*(1-b*T(t))-c1*N(t)*T(t)-D(t)*T(t)-KT*M(t)*T(t), #tumor cells
    L'(t) = -m*L(t) - q*L(t)*T(t)-ucte*L(t)^2+r2*C(t)*T(t) + pI*L(t)*I(t)/(gI+I(t)) + u1(t) - KL*M(t)*L(t), # tumor-specific effector cells, T-celss
    N'(t) = alpha1 - f*N(t)+ g*T(t)*N(t)/(h + T(t)) - p*N(t)*T(t) - KN*M(t)*N(t), # non-specific effector cells, NK cells
    C'(t) = alpha2 - beta*C(t) - KC*M(t)*C(t), #circulating lymphocytes
    I'(t) = pt*T(t)*L(t)/(gt + T(t)) +  w*L(t)*I(t) - muI*I(t) + u2(t), # IL-2, VI = u2 aplicación directa, terapia de IL2
    M'(t) = - gamma*M(t) + u1(t), #chemotherapy drug, terapia/aplicación de quimio, u1 = VM
    y1(t) = L(t) + N(t) + C(t),
    y2(t) = M(t),
    y3(t) = I(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


#
using SIAN, Logging

ode = @ODEmodel(
    T'(t) = a*T(t)*(1-b*T(t))-c1*N(t)*T(t)-D(t)*T(t)-KT*M(t)*T(t), #tumor cells
    L'(t) = -m*L(t) - q*L(t)*T(t)-ucte*L(t)^2+r2*C(t)*T(t) + pI*L(t)*I(t)/(gI+I(t)) + u1(t) - KL*M(t)*L(t), # tumor-specific effector cells, T-celss
    N'(t) = alpha1 - f*N(t)+ g*T(t)*N(t)/(h + T(t)) - p*N(t)*T(t) - KN*M(t)*N(t), # non-specific effector cells, NK cells
    C'(t) = alpha2 - beta*C(t) - KC*M(t)*C(t), #circulating lymphocytes
    I'(t) = pt*T(t)*L(t)/(gt + T(t)) +  w*L(t)*I(t) - muI*I(t) + u2(t), # IL-2, VI = u2 aplicación directa, terapia de IL2
    M'(t) = - gamma*M(t) + u1(t), #chemotherapy drug, terapia/aplicación de quimio, u1 = VM
    y1(t) = L(t) + N(t),
    Y4(t) = C(t),
    y2(t) = M(t),
    y3(t) = I(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
