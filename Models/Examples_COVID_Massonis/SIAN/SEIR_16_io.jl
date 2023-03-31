using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -beta*I(t)*S(t),
    E'(t) = beta*I(t)*S(t)-epsilon*E(t),
    I'(t) = epsilon*E(t)-(rho+mu)*I(t),
    R'(t) = rho*I(t)-d*R(t),
    y1(t) = mu*I(t)
)

println(find_ioequations(ode))
    
using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -beta*I(t)*S(t),
    E'(t) = beta*I(t)*S(t)-epsilon*E(t),
    I'(t) = epsilon*E(t)-(rho+mu)*I(t),
    R'(t) = rho*I(t) - R(t),
    y1(t) = mu*I(t)
)

println(find_ioequations(ode))
