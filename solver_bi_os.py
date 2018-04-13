from fenics import *
import numpy as np
import time
import os
import matplotlib.pyplot as plt

M_i  = 3
M_e  = 3
sigma = 1
c_1 = 200
a_1 = 0.1
c_2 = 200
c_3 = 1
b = 1
num_jacobi_iter = 1

Nx = 120
Ny = Nx

degree = 2
two_d = 1
T = 0.075
num_steps = 400

plotting = True
advance_method = "OS"           #"FE" if ForwardEuler and "OS" if OperatorSplitting
dt = T/float(num_steps)
print dt
start = time.time()

class ODEsolver:
    def __init__(self, T, num_steps, plotting, num_elements, advance_method):
        self.T, self.num_steps = T, num_steps
        self.plotting = plotting
        self.v_1_vector = v_0.vector()
        self.w_1_vector = w_0.vector()
        self.u_1_vector = u_0.vector()
        self.dt = T/float(num_steps)
        self.num_elements = num_elements
        self.tid = 0
        self.advance_method = advance_method

    def solver(self):
        filename = "bi_%sd_%s_%s_%s" % (2**two_d, self.num_elements, \
                                self.num_steps, self.advance_method)
        # time stepping
        for n in range(0, self.num_steps ):
            self.n = n
            #if n > 6:
            #    break
            if n == 0 and advance_method == "FE":
                # loser for u_e i forste tidssteg med conjugate gradient metoden
                solve(A2, u.vector(), -1*A1*self.v_1_vector, "cg")
                self.u_1_vector = u.vector()
                print self.u_1_vector.get_local()
            v.vector()[:], w.vector()[:], u.vector()[:] = self.advance()
            self.w_1_vector = w.vector()
            self.v_1_vector = v.vector()
            self.u_1_vector = u.vector()
            self.tid += self.dt
            if n%20 == 0:
                  print self.tid
            # break if numerical solution diverges
            if abs(np.sum(self.v_1_vector.get_local())/len(self.v_1_vector.get_local())) > 10:
                print "break"
                break
            if self.plotting and n%int(num_steps/4) == 0:
                plot(v)
                plt.show()
        output_file = XDMFFile(filename + ".xdmf")
        output_file.write_checkpoint(v,'v')
        File(filename + '_mesh.xml.gz') << v.function_space().mesh()

class OperatorSplitting(ODEsolver):
    """ Bruker Godunov-splitting og baklengs Euler
    """
    def advance(self):
        v_1_vector = self.v_1_vector
        w_1_vector = self.w_1_vector
        dt = self.dt
        # step 1
        v_n_s1 = v_1_vector + dt*(c_1*v_1_vector*\
                (v_1_vector - a_1)*(1 - v_1_vector) - c_2*v_1_vector*w_1_vector)

        w_n_s1 = dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector
        # step 2
        # define variational problem for implicit scheme
        _v.vector()[:] = v_n_s1
        v_var, u_var = split(u_v)
        F = (v_var - _v)/dt*psi1*dx + inner(M_i*grad(v_var),grad(psi1))*dx + inner(M_i*grad(u_var),grad(psi1))*dx + \
            inner(M_i*grad(v_var), grad(psi2))*dx + inner((M_i + M_e)*grad(u_var), grad(psi2))*dx
        solve(F == 0, u_v)
        u_v_array = u_v.vector().get_local()
        v_split, u_split = u_v.split(deepcopy=True)
        return v_split.vector(), w_n_s1, u_split.vector()

vtk_v = File("v_bidomain.pvd")

# create Expressions
if two_d:
    I_v_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = 2)
    I_w_expression = Expression("x[0] < 0.2 && x[1] < 0.2 ? 0: 0", degree = 1)
    I_u_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = 2)
    mesh_mix = UnitSquareMesh(Nx, Ny)

else:
    I_v_expression = Expression("x[0] < 0.2 ? 1: 0, 0", degree = 1)
    I_w_expression = Expression("x[0] < 0.2 ? 0: 0", degree = 2)
    I_u_expression = Expression("x[0] < 0.2 ? 1: 0", degree = 2)
    mesh_mix = UnitIntervalMesh(Nx)

#create function space
V_fe_mix = FiniteElement("CG", mesh_mix.ufl_cell(), degree)
V = FunctionSpace(mesh_mix, MixedElement(V_fe_mix, V_fe_mix))
W = FunctionSpace(mesh_mix, V_fe_mix)
w = Function(W)
u = Function(W)
v = Function(W)
_v = Function(W)
_u = Function(W)
u_v = Function(V)

(psi1, psi2) = TestFunctions(V)
#interpolate initial condition
v_0 = interpolate(I_v_expression, W)
w_0 = interpolate(I_w_expression, W)
u_0 = interpolate(I_u_expression, W)

object_os = OperatorSplitting(T = T, num_steps = num_steps, plotting = plotting,
                                num_elements = Nx, advance_method = advance_method)

solution_os = object_os.solver()

end = time.time()
print end - start
#interactive()
