from fenics import *
import numpy as np
import time
import solver_monodomene as solver_mono
import solver_bidomene as solver_bi
from math import log

def compute_rates(dt_values, E_values):
    m = len(dt_values)
    r = [log(E_values[i-1]/E_values[i])/log(dt_values[i-1]/dt_values[i]) \
            for i in range(1, m, 1)]
    return r

def table(num_steps, num_elements, advance_method, two_d, mono, T = 0.1, plotting = False,
        sigma = 0.1, M_i  = 0.1, M_e = 0.1):
    degree = 1
    degree_ref = 1
    dt_values = [T/float(timestep) for timestep in num_steps]
    if mono:
        if two_d:
            filename_ref = "2d_600_16000_False_True_FE_0.100"
        else:
            filename_ref = "1d_10000_2100000_False_True_FE_0.100"
    else:
        if two_d:
            filename_ref = "2d_600_16000_False_True_FE_0.100"
        else:
            filename_ref = "1d_10000_2100000_False_True_FE_0.100"
    Error_values = []
    vectors_v = []
    for j in range(len(num_steps)):
        if mono:
            if advance_method == "OS":
                object_ = solver_mono.OperatorSplitting(T = T, num_steps = num_steps[j], plotting = plotting,
                                    num_elements = num_elements[j], system = True, heat = False,
                                     advance_method = advance_method, two_d = two_d,
                                     sigma = sigma, degree = degree)
            else:
                object_ = solver_mono.ForwardEuler(T = T, num_steps = num_steps[j], plotting = plotting,
                                    num_elements = num_elements[j], system = True, heat = False,
                                     advance_method = advance_method, two_d = two_d,
                                     sigma = sigma, degree = degree)
        else:
            "hei"
        object_.solver()
        vectors_v.append(object_.v)
    input_file = XDMFFile(filename_ref +'.xdmf')
    mesh_v = Mesh(filename_ref +'_mesh.xml.gz')
    V_ref_v = FunctionSpace(mesh_v, 'P', degree_ref)
    v_ref = Function(V_ref_v)
    input_file.read_checkpoint(v_ref, "v")
    vectors_v.append(v_ref)
    for j in range(0, len(vectors_v) - 1, 1):
        v = vectors_v[j]
        v_ref_interpolate = interpolate(v_ref, v.ufl_function_space())
        error_L2 = errornorm(v_ref_interpolate, v, "L2")
        Error_values.append(error_L2)
    conv_rates = compute_rates(dt_values, Error_values)
    print '    dx            dt          L2 error  conv. rate'
    print "%6.6f% 6.6f% 10.4f" % (1/float(num_elements[0]), dt_values[0],
                                Error_values[0])
    for j in range(len(conv_rates)):
        print "%6.6f% 6.6f% 10.4f%12.4f" % (1/float(num_elements[j+1]), dt_values[j+1],
                                Error_values[j+1],conv_rates[j])

if __name__ = "main":
    table(num_steps = [500], num_elements = [25], advance_method = "FE", two_d = False,
    mono = True, plotting = True)
