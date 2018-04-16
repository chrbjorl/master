from fenics import *
import numpy as np
import time
from solver_monodomene import *
from math import log

num_experiments = 1

sigma = 0.03
T = 0.2

model = "mono"
plotting = False
num_steps_ref = 10000
num_steps = [400 for i in range(11)]
#num_steps = [50*2**i for i in range(8)]

num_elements_ref = 500
num_elements =[40 + i for i in range(11)]
#num_elements =[5*2**i for i in range(8)]
method_ref = "FE"
method = "FE"

num_elements.append(num_elements_ref)
num_steps.append(num_steps_ref)


heat = False
system = False
two_d = 0
degree = 1
degree_ref = 1

dt_values = [T/float(timestep) for timestep in num_steps[0:-1]]
Error_values = []


def compute_rates(dt_values, E_values):
    m = len(dt_values)
    r = [log(E_values[i-1]/E_values[i])/log(dt_values[i-1]/dt_values[i]) \
            for i in range(1, m, 1)]
    return r

vectors = []
for j in range(len(num_steps)):
    if j == len(num_steps) - 1:
        method = method_ref
    if model == "bi":
        filename = "bi_%sd_%s_%s_%s_%s" % (2**two_d, num_elements[j], num_steps[j],
                                        method, num_jacobi_iter[j])
        if method == "FE":
            object_fe = ForwardEuler(T = T, num_steps = num_steps[j], plotting = plotting,\
                                        num_elements = num_elements[j], system = system, heat = heat,
                                        advance_method = advance_method)
            solution_fe = object_fe.solver()

        else:
            object_os = OperatorSplitting(T = T, num_steps = num_steps[j], plotting = plotting,\
                                        num_elements = num_elements[j], advance_method = advance_method)
            solution_os = object_os.solver()

    if model == "mono":
        filename = "%sd_%s_%s_%s_%s_%s_%1.3f" % (2**two_d, num_elements[j], num_steps[j],
                                                        heat, system, method, sigma)
        if method == "FE":
            object_fe = ForwardEuler(T = T, num_steps = num_steps[j], plotting = plotting,
                                        num_elements = num_elements[j], system = system, heat = heat,
                                        advance_method = advance_method, two_d = two_d, sigma = sigma)
            solution_fe = object_fe.solver()
        else:
            object_fe = ForwardEuler(T = T, num_steps = num_steps[j], plotting = plotting,
                                        num_elements = num_elements[j], system = system, heat = heat,
                                        advance_method = advance_method, two_d = two_d, sigma = sigma)
            solution_os = object_os.solver()
    input_file = XDMFFile(filename+'.xdmf')
    print filename
    mesh = Mesh(filename+'_mesh.xml.gz')
    V_ref = FunctionSpace(mesh, 'P', degree)
    u = Function(V_ref)
    input_file.read_checkpoint(u,'v')
    vectors.append(u)

v_ref = vectors[-1]
for j in range(0, len(vectors) - 1, 1):
    v = vectors[j]
    v_ref_interpolate = interpolate(v_ref, v.ufl_function_space())
    #print v_ref_interpolate.vector().get_local()
    error_L2 = errornorm(v_ref_interpolate, v, "L2")
    Error_values.append(error_L2)

print Error_values

conv_rates = compute_rates(dt_values, Error_values)
print '    dt       dx    L2 error  conv. rate'
for j in range(len(conv_rates)):
    print "%6.6f%7.8f%10.6f%12.4f" % (dt_values[j+1], num_elements[j+1],
                                Error_values[j+1],conv_rates[j])
