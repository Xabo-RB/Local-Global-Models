#Saccomani, M. P., & Thomaseth, K. (2019). Calculating all multiple parameter solutions of ODE models to avoid biological 
#misinterpretations. Mathematical Biosciences and Engineering, 16(6), 6438-6453.

using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -(k01 + k21 + k31)*x1(t) + k12*x2(t) + u(t),
    x2'(t) = k21*x1(t) - (k12 + k32)*x2(t) + k23*x3(t),
    I'(t) = k31*x1(t) + k32*x2(t) - (k23+k03)*x3(t),
    y1(t) = x2(t),
    y2(t) = x3(t)
)

println(find_ioequations(ode))
