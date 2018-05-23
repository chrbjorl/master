import solver_bidomene as solver_bi
from efficiency import error_
from fenics import *


Nt_os = 70
Nx_os = 40
num_jacobi_iter = [1, 1, 3, 1]
jacobi_ = [False, True, True, True]
Nt_fe = [70, 70, 70, 250]
Nx_fe = [36 for i in range(4)]
two_d = True
mono = False

object_os = solver_bi.OperatorSplitting(T = 0.1, num_steps = Nt_os,
                        plotting_v = False, num_elements = Nx_os,  advance_method = "OS",
                         two_d = two_d)
object_os.solver()
time_os = object_os.end - object_os.start
error_os = error_(two_d = two_d, object_ = object_os, mono = mono)

errors_fe = []
time_fe = []
for j in range(len(num_jacobi_iter)):
    object_fe = solver_bi.ForwardEuler(T = 0.1, num_steps = Nt_fe[j],
                        plotting_v = False, num_elements = Nx_fe[j],
                         advance_method = "FE", two_d = two_d,
                         num_jacobi_iter = num_jacobi_iter[j], jacobi_ = jacobi_[j])
    object_fe.solver()
    errors_fe.append(error_(two_d = two_d, object_ = object_fe, mono = mono))
    time_fe.append(object_fe.end - object_fe.start)

print errors_fe
print time_fe


print "Method for solving elliptic part            CPU-time     dx      dt     error"
print "conjugate gradient method until convergence   %1.2f  %1.6f %1.6f %1.4f" %\
        (time_fe[0], 1/float(Nx_fe[0]), 1/float(Nt_fe[0]), errors_fe[0])
print "one jacobi iteration                          %1.2f  %1.6f %1.6f %1.4f" %\
     (time_fe[1], 1/float(Nx_fe[1]), 1/(float(Nt_fe[1])*10), errors_fe[1])
print "three jacobi iteration                        %1.2f  %1.6f %1.6f %1.4f" %\
     (time_fe[2], 1/float(Nx_fe[2]), 1/(float(Nt_fe[2])*10), errors_fe[2])
print "one jacobi iteration                          %1.2f  %1.6f %1.6f %1.4f" %\
         (time_fe[3], 1/float(Nx_fe[3]), 1/(float(Nt_fe[3])*10), errors_fe[3])
print "operator splitting                            %1.2f  %1.6f %1.6f %1.4f" %\
         (time_os, 1/float(Nx_os), 1/(float(10*Nt_os)*10), error_os)
