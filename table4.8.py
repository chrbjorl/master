import solver_monodomene as solver_mono
from efficiency import error_
from fenics import *

Nt_os = 36
Nx_os = 35
Nt_fe = 180
Nx_fe = 44
two_d = True
mono = True

object_os = solver_mono.OperatorSplitting(T = 0.1, num_steps = Nt_os,
                        plotting = False, num_elements = Nx_os, system = True,
                        heat = False, advance_method = "OS", two_d = two_d)
object_os.solver()
time_os = object_os.end - object_os.start
error_os = error_(object_ = object_os, mono = mono, two_d = two_d)

object_fe = solver_mono.ForwardEuler(T = 0.1, num_steps = Nt_fe,
                        plotting = False, num_elements = Nx_fe, system = True,
                        heat = False, advance_method = "FE", two_d = two_d)
object_fe.solver()
time_fe = object_fe.end - object_fe.start
error_fe = error_(object_ = object_fe, mono = mono, two_d = two_d)

print "numerical scheme   CPU-time   dx       dt      error"
print "operator splitting     %1.2f %1.6f %1.6f %1.4f" % (time_os,
                    1/float(Nx_os), 1/float(10*Nt_os), error_os*100)
print "explicit               %1.2f %1.6f %1.6f %1.4f" % (time_fe,
                    1/float(Nx_fe), 1/float(10*Nt_fe), error_fe*100)
