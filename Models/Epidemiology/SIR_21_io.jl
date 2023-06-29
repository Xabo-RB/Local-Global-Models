#Roda, W. C., Varughese, M. B., Han, D., & Li, M. Y. (2020). Why is it difficult to accurately predict the COVID-19 epidemic?.
# Infectious disease modelling, 5, 271-281.


using StructuralIdentifiability


ode = @ODEmodel(
    N'(t) = 0,
    S'(t) = -beta*S(t)*I(t)/N(t)-pp*S(t)+q*C(t),
    I'(t) = beta*S(t)*I(t)/N(t)-(r+mu)*I(t),
    R'(t) = r*I(t),
    C'(t) = pp*S(t)-q*C(t),
    D'(t) = mu*I(t),
    y1(t) = N(t),
    y2(t) = D(t),
    y3(t) = C(t)
)

#@time println(assess_identifiability(ode))

println(find_ioequations(ode))
