using StructuralIdentifiability

ode = @ODEmodel(
    S'(t) = -lambda * phi * (S(t)-T(t)) + s*(1-S(t)),
    T'(t) =  phi * (S(t) - T(t)) + s*(1-T(t)) - k * (T(t)^5)*(L^5),
    A'(t) =  k * (T(t)^5)*(L^5) - ki*A(t),
    y1(t) = S(t)/(lambda+1) + ((T(t)+A(t))*lambda/(lambda+1))
)

println(find_ioequations(ode))
