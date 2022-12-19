#Hu, X., Ke, G., & Jang, S. R. J. (2019). 
#Modeling pancreatic cancer dynamics with immunotherapy. Bulletin of mathematical biology, 81(6), 1885-1915.


using SIAN, Logging

ode = @ODEmodel(
    S'(t) = -lambda * phi * (S(t)-T(t)) + s*(1-S(t)),
    T'(t) =  phi * (S(t) - T(t)) + s*(1-T(t)) - k * (T(t)^5)*(L^5),
    A'(t) =  k * (T(t)^5)*(L^5) - ki*A(t),
    y1(t) = S(t)/(lambda+1) + ((T(t)+A(t))*lambda/(lambda+1))
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
