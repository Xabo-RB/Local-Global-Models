using StructuralIdentifiability

ode = @ODEmodel(
    x12'(t) = (1-Theta1-Theta2)*b - (m1*Lambda1 + m2*Lambda2 + mu)*x12 + (nu1+tau)*x22(t) + (nu2+tau)*x3(t) + tau*x4(t),
    x22'(t) = Theta1*b + m1*Lambda1*x12 + nu2*x4(t) -((1-Pi2)*m2*Lambda2+nu1+mu+c1+tau)*x22(t),
    x3'(t) = Theta2*b + m2*Lambda2*x12 + nu1*x4(t) -((1-Pi1)*m1*Lambda1+nu2+mu+tau)*x3(t),
    x4'(t) = (1-Pi1)*m1*Lambda1*x3(t) + (1-Pi2)*m2*Lambda2*x22(t) - (nu1+nu2+mu+c1+tau)*x4(t),
    y1(t) = x12+x22(t)+x3(t)+x4(t),
    y2(t) = x22(t)+x4(t),
    y3(t) = x3(t)+x4(t)
)

@time println(assess_identifiability(ode))