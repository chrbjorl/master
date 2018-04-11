from fenics import *
import numpy as np
import time
from math import log

num_experiments = 1

sigma = 1


model = "bi"

num_steps_ref = 2000
num_steps = [2000]

num_elements_ref = 30
num_elements = [30]

method_ref = "FE"
method = "OS"

num_elements.append(num_elements_ref)
num_steps.append(num_steps_ref)


heat = False
two_d = False
degree = 1
degree_ref = 1

error_rates = len(num_elements) > 2.1
T = 1


def compute_rates(dt_values, E_values):
    m = len(dt_values)
    r = [log(E_values[i-1]/E_values[i])/log(dt_values[i-1]/dt_values[i]) \
            for i in range(1, m, 1)]
    return r

def read(filename, num_elements, degree):
    infile = open(filename, "r")
    num_points = (2*num_elements + 1)**2 if degree == 2 else num_elements + 1
    u = np.zeros(num_points)
    for line in infile:
        words = line.split()
        if words[0] != "numerical" and words[0] != "exact":
            for w in range(len(words)):
                u[w] = float(words[w])
    return u


dt_values = [T/float(timestep) for timestep in num_steps[0:-1]]
Error_values = []
l2_times = 0.5

u_elements = []
for j in range(len(num_steps[0:-1])):
    if model == "bi":
        filename = "bi_%sd_%s_%s_%s_%s_%1.3f.dat" % (degree, num_elements[j], num_steps[j],
                                        heat, method, sigma)
    else:
        filename = "%sd_%s_%s_%s_%s_%1.3f.dat" % (degree, num_elements[j], num_steps[j],
                                        heat, method, sigma)
    u = read(filename, num_elements[j], degree)
    u_elements.append(u)
    print filename

if model == "bi":
    filename = "bi_%sd_%s_%s_%s_%s_%1.3f.dat" % (degree_ref, num_elements_ref, num_steps_ref,
                                    heat, method_ref, sigma)
else:
    filename = "%sd_%s_%s_%s_%s_%1.3f.dat" % (degree_ref, num_elements_ref, num_steps_ref,
                                        heat, method_ref, sigma)
print filename
u_elements.append(read(filename, num_elements_ref, degree_ref))
if degree_ref == 2:
    mesh_ref = UnitSquareMesh(num_elements_ref, num_elements_ref)
if degree_ref == 1:
    mesh_ref = UnitIntervalMesh(num_elements_ref)

V_ref = FunctionSpace(mesh_ref, "P", degree_ref)
v_ref = Function(V_ref)
v_ref.vector()[:] = u_elements[-1][:]

def l2_error(a,b, length):
    return np.sqrt(sum((a-b)**2)/length)

for j in range(0, len(u_elements) - 1, 1):
    if degree == 2:
        mesh = UnitSquareMesh(num_elements[j], num_elements[j])
    if degree == 1:
        mesh = UnitIntervalMesh(num_elements[j])
    V = FunctionSpace(mesh, "P", degree)
    v = Function(V)
    v.vector()[:] = u_elements[j][:]
    v_ref = interpolate(v_ref, v.ufl_function_space())
    error_L2 = errornorm(v, v_ref, "L2")
    print error_L2


Error_values = [0.2839116197, 0.2314466729, 0.0193766096, 0.0061085585, 0.0016177889]

print compute_rates(dt_values, Error_values)
