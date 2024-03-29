using StructuralIdentifiability

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t),
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

println(find_ioequations(ode))


using StructuralIdentifiability

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t),
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

println(find_ioequations(ode))


using StructuralIdentifiability

ode = @ODEmodel(
    x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t),
    x2'(t) = x1(t)-k12*x2(t), #elimino k21
    x3'(t) = k31*x1(t)-k13*x3(t),
    x4'(t) = k41*x1(t)-k14*x4(t),
    y1(t) = x1(t)
)

println(find_ioequations(ode))
