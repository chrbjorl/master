from fenics import *
import numpy as np
import time
from solver_monodomene import *
from math import log

sigma = 0.1
T = 0.1

system = True
two_d = 1
degree = 1
degree_ref = 1
plotting = False
heat = False
method_ref = "FE"
advance_method = "OS"
if two_d:
    num_elements_ref = 600
    num_steps_ref = 16000

else:
    num_elements_ref = 10000
    num_steps_ref = 2100000

num_steps = [180]#*4**i for i in range(3)]
num_elements = [60*2**i for i in range(3)]

num_elements.append(num_elements_ref)
num_steps.append(num_steps_ref)

dt_values = [T/float(timestep) for timestep in num_steps[0:-1]]
Error_values = []

def compute_rates(dt_values, E_values):
    m = len(dt_values)
    r = [log(E_values[i-1]/E_values[i])/log(dt_values[i-1]/dt_values[i]) \
            for i in range(1, m, 1)]
    return r

vectors_v = []
for j in range(len(num_steps)):
    if abs(j - 1) < 1e-08:
        filename1 = filename
    if abs(j -(len(num_steps) - 1)) < 0.1:
        advance_method = method_ref
    filename = "%sd_%s_%s_%s_%s_%s_%1.3f" % (2**two_d, num_elements[j], num_steps[j],
                                                        heat, system, advance_method, sigma)

    if abs(j -(len(num_steps) - 1)) > 0.1:
        print filename
        if advance_method == "OS":
            object_os = OperatorSplitting(T = T, num_steps = num_steps[j], plotting = plotting,
                                    num_elements = num_elements[j], system = system, heat = heat,
                                     advance_method = advance_method, two_d = two_d,
                                     sigma = sigma, degree = degree)
            object_os.solver()
        else:
            object_fe = ForwardEuler(T = T, num_steps = num_steps[j], plotting = plotting,
                                    num_elements = num_elements[j], system = system, heat = heat,
                                     advance_method = advance_method, two_d = two_d,
                                     sigma = sigma, degree = degree)
            object_fe.solver()
    input_file = XDMFFile(filename+'.xdmf')
    mesh_v = Mesh(filename +'_mesh.xml.gz')
    V_ref_v = FunctionSpace(mesh_v, 'P', degree)
    v = Function(V_ref_v)
    input_file.read_checkpoint(v, "v")
    vectors_v.append(v)

v_ref = vectors_v[-1]
for j in range(0, len(vectors_v) - 1, 1):
    v = vectors_v[j]
    v_ref_interpolate = interpolate(v_ref, v.ufl_function_space())
    error_L2 = errornorm(v_ref_interpolate, v, "L2")
    Error_values.append(error_L2)

#print Error_values_u
print Error_values
print filename
conv_rates = compute_rates(dt_values, Error_values)
print '    dx            dt          L2 error  conv. rate'
print "%6.6f% 6.6f% 10.4f" % (1/float(num_elements[0]), dt_values[0],
                                Error_values[0])
for j in range(len(conv_rates)):
    print "%6.6f% 6.6f% 10.4f%12.4f" % (1/float(num_elements[j+1]), dt_values[j+1],
                                Error_values[j+1],conv_rates[j])

print filename1
