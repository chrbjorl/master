from fenics import *

def error_(two_d, object_, mono, degree = 1):
    if mono:
        if two_d:
            filename_ref = "2d_600_16000_False_True_FE_0.100"
        else:
            filename_ref = "1d_10000_2100000_False_True_FE_0.100"
    else:
        if two_d:
            filename_ref = "bi_2d_500_500_OS"
        else:
            filename_ref = "bi_1d_10000_40000_OS"
    input_file = XDMFFile(filename_ref +'.xdmf')
    if mono:
        mesh_v = Mesh(filename_ref +'_mesh.xml.gz')
    else:
        mesh_v = Mesh(filename_ref +'_mesh_v.xml.gz')
    V_ref_v = FunctionSpace(mesh_v, 'P', degree)
    v_ref = Function(V_ref_v)
    input_file.read_checkpoint(v_ref, "v")
    v_ref_interpolate = interpolate(v_ref, object_.v.ufl_function_space())
    error_L2 = errornorm(v_ref_interpolate, object_.v, "L2")
    return error_L2
