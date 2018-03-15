from fenics import *
import numpy as np
import time
import numpy.linalg as LA
from scipy.linalg import solve_banded
from math import log

num_experiments = 1

num_steps_ref = 500
num_steps = [500]

num_elements_ref = 10
num_elements = [2]

method_ref = "OS"
method = "OS"

num_elements.append(num_elements_ref)
num_steps.append(num_steps_ref)

heat = True
degree = 2
error_rates = len(num_elements) > 2.1
T = 1
l2_time = 0.5

def compute_rates(dt_values, E_values):
    m = len(dt_values)
    r = [log(E_values[i-1]/E_values[i])/log(dt_values[i-1]/dt_values[i]) \
            for i in range(1, m, 1)]
    return r

def read(filename, num_elements):
    infile = open(filename, "r")
    u = np.zeros((2*num_elements + 1)**2)
    for line in infile:
        words = line.split()
        if words[0] != "numerical" and words[0] != "exact":
            for w in range(len(words)):
                u[w] = float(words[w])
    return u


#trenger ikke a ha med referanselosning i dt_values
dt_values = [T/float(timestep) for timestep in num_steps[0:-1]]
Error_values = []
l2_times = 0.5

u_elements = []
for num_element, timestep in zip(num_elements[0:-1], num_steps[0:-1]):
    filename = "%sd_%s_%s_%1.9f_%s_%s.dat" % (degree, num_element,
                                                    timestep, l2_time, heat, method)
    u = read(filename, num_elements_ref)
    u_elements.append(u)

filename = "%sd_%s_%s_%1.9f_%s_%s.dat" % (degree, num_elements_ref, num_steps_ref,
                                                l2_time, heat, method_ref)

u_elements.append(read(filename, num_elements_ref))
mesh_ref = UnitSquareMesh(num_elements_ref, num_elements_ref)
V_ref = FunctionSpace(mesh_ref, "P", degree)
v_ref = Function(V_ref)
v_ref.vector()[:] = u_elements[-1][:]

print u_elements[0]


for j in range(0, len(u_elements) - 1, 1):
    mesh = UnitSquareMesh(num_elements[j], num_elements[j])
    V = FunctionSpace(mesh, "P", degree)
    v = Function(V)

    v.vector()[:] = u_elements[j][:]
    v_ref = interpolate(v_ref, v.ufl_function_space())
    #plot(v_ref)
    error_L2 = errornorm(v, v_ref, "L2")
    Error_values.append(error_L2)
    #s_error += "%10.9f" %

print Error_values
"""
outfile = open("errors.dat", "a")
outfile.write("\n")
s = "#timesteps1 = %s, #timesteps2 = %s \n" % (timesteps[0], timesteps[1])
outfile.write(s)
s = "#elements1 = %s, #elements2 = %s \n" % (num_elements[0], num_elements[1])
outfile.write(s)
s = "#method1 = %s, #tmethod2 = %s \n" % (method, method_ref)
outfile.write(s)
s = "error = %10.8f \n" % (error_L2)
outfile.write(s)
if error_rates:
    rates = compute_rates(dt_values, Error_values)
    print rates
"""
