using SIAN, Logging
ode = @ODEmodel(
    beta'(t) = -kbeta*beta(t),
    cry'(t) =  -kcryOH*cry(t)-kcrybeta*cry(t),
    zea'(t) = -kzea*zea(t),
    beta10'(t) = kbeta*beta(t)+kcryOH*cry(t)-kbeta10*beta10(t),
    OHbeta10'(t) = kcrybeta*cry(t) + kzea*zea(t) - kOHbeta10*OHbeta10(t),
    betaio'(t) = kbeta*beta(t) + kcrybeta*cry(t) + kbeta10*beta10(t),
    OHbetaio'(t) = kcryOH*cry(t) + kzea*zea(t) + kOHbeta10*OHbeta10(t),
    y1(t) = beta(t),
    y2(t) = beta10(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))


