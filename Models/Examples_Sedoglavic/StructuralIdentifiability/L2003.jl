using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -theta1*(x1(t)/theta2-x2(t)/theta3)/(1+x1(t)/theta2+x2(t)/theta3)+ theta4*(u(t)-x1(t))/(1+u(t)/theta5+x1(t)/theta5+x1(t)*u(t)/theta5^2),
    x2'(t) = theta1*(x1(t)/theta2-x2(t)/theta3)/(1+x1(t)/theta2+x2(t)/theta3)-theta6*(x2(t)/theta7-theta8)/(1+x2(t)/theta7+theta8),
    y1(t) = x2(t)
)

@time println(assess_identifiability(ode))