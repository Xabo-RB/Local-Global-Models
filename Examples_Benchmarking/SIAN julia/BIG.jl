# Diane T Finegood, Luisa Scaglia, and Susan Bonner-Weir. Dynamics of β-cell mass in the growing rat
# pancreas: estimation with a simple mathematical model. Diabetes, 44(3):249–256, 1995

using Logging, SIAN



ode = @ODEmodel(
    G'(t) = u1(t)-(c+si*II(t))*G(t),
    beta'(t) = beta(t)* ( ( (1/(48*24*60)) /(1+((8/G(t))^2)) ) - ( (1/(40*24*60))/(1+((G(t)/5)^8)))),
    II'(t) = p1*beta(t)*(G(t)^2/(alpha^2+G(t)^2)) - gm*II(t),
    y1(t) = G(t)
)


@time res = identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 0, nthrds = 1)

println(res)









