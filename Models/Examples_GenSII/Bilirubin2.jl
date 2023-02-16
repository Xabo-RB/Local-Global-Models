#Bellu, G., Saccomani, M. P., Audoly, S., & D’Angiò, L. (2007). DAISY: A new software tool to test global identifiability of biological and physiological systems. 
#Computer methods and programs in biomedicine, 88(1), 52-61.

using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t),
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#PRUEBAS
# todo NI incluido k01
using SIAN, Logging
ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t), #elimino U
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k21 NI, k12 SGI
using SIAN, Logging
ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t),
  x2'(t) = x1(t)-k12*x2(t), #elimino k21
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k12 y k21 a NI y k01 a NI
using SIAN, Logging
ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t),
  x2'(t) = k21*x1(t)-x2(t), #elimino k12
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k01 SLI
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+x2(t)+k13*x3(t)+k14*x4(t)+u(t), #elimino k12
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k01 SLI
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t), #elimino k21
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k01 SLI
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+x2(t)+k13*x3(t)+x4(t)+u(t), #eliminé k12 y k14
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k41 y k14 y x4 a SGI/ k21 y k31 y k01 a NI
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t),
  x2'(t) = x1(t)-k12*x2(t), #eliminé k21
  x3'(t) = x1(t)-k13*x3(t), #eliminé k31
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t)
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k12 y k21 y x2 a SGI
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t)+u(t),
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t),
  y2(t) = x2(t) #como en bilirubin 1
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))

#k12 y k21 y x2 a SGI, k01 y todo SLI a SGI
using SIAN, Logging

ode = @ODEmodel(
  x1'(t) = -(k21+k31+k41+k01)*x1(t)+k12*x2(t)+k13*x3(t)+k14*x4(t),
  x2'(t) = k21*x1(t)-k12*x2(t),
  x3'(t) = k31*x1(t)-k13*x3(t),
  x4'(t) = k41*x1(t)-k14*x4(t),
  y1(t) = x1(t),
  y2(t) = x2(t) #como en bilirubin 1
)

@time println(identifiability_ode(ode, get_parameters(ode); p = 0.99, p_mod = 2^29 - 3, infolevel = 10, nthrds = 1))
