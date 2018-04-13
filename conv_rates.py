from fenics import *
import numpy as np
import time
from math import log

num_experiments = 1

sigma = 1
T = 1

model = "bi"

num_steps_ref = 200
num_steps = [50, 100]

num_elements_ref = 40
num_elements =[10, 20]
method_ref = "OS"
method = "OS"

num_elements.append(num_elements_ref)
num_steps.append(num_steps_ref)


heat = False
two_d = 1
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
    if model == "bi":
        filename = "bi_%sd_%s_%s_%s" % (2**two_d, num_elements[j], num_steps[j],
                                        method)
    else:
        filename = "%sd_%s_%s_%s_%s_%1.3f" % (2**two_d, num_elements[j], num_steps[j],
                                        heat, method, sigma)
    input_file = XDMFFile(filename+'.xdmf')
    mesh = Mesh(filename+'_mesh.xml.gz')
    V_ref = FunctionSpace(mesh, 'P', degree)
    u = Function(V_ref)
    input_file.read_checkpoint(u,'v')
    vectors.append(u)

for j in range(0, len(vectors) - 1, 1):
    v_ref = vectors[-1]
    v = interpolate(vectors[j], v_ref.ufl_function_space())
    error_L2 = errornorm(v, v_ref, "L2")
    Error_values.append(error_L2)
    array = v.vector().get_local()
    print error_L2


print Error_values
print compute_rates(dt_values, Error_values)
