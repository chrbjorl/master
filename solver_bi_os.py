from fenics import *
import numpy as np
import time
import scipy.sparse
import os
#os.sys.path.append("./cbc.block/block")
#from block import *
#from block.iterative import *
#from block.algebraic.trilinos import ML


M_i  = 3
M_e  = 3
sigma = 1
c_1 = 200
a_1 = 0.1
c_2 = 200
c_3 = 1
b = 1
num_jacobi_iter = 1

Nx = 30
Ny = Nx

degree = 1
two_d = 0
T = 0.075
num_steps = 2000

system_ = True            # fitzhugh-nagumo if system = True
heat = False
plotting = False
advance_method = "OS"           #"FE" if ForwardEuler and "OS" if OperatorSplitting
dt = T/float(num_steps)
print dt
def jacobi(u_0, K, v_fe):
    u_old = u_0
    product1 = A1*v_fe
    for k in range(K):
        #product1 = -1*A1*v_fe
        #u_new = -0.01041667*(product1 + R*u_old)
        u_new = u_old + alfa_jacobi*(product1 - A2*u_old)
        u_old = u_new
    return u_new
start = time.time()

class ODEsolver:
    def __init__(self, T, num_steps, plotting, ref, num_elements,
                system_, heat, exact_expression, advance_method):
        self.T, self.num_steps = T, num_steps
        self.plotting, self.ref = plotting, ref
        self.v_1_vector = v_0.vector()
        self.w_1_vector = w_0.vector()
        self.u_1_vector = u_0.vector()
        self.dt = T/float(num_steps)
        self.num_elements = num_elements
        self.tid = 0
        self.system_ = system_
        self.heat = heat
        self.exact_expression = exact_expression
        self.advance_method = advance_method


    def solver(self):
        ref = self.ref
        filename_referanse = "bi_%sd_%s_%s_%s_%s_%1.3f.dat" % (2**two_d, self.num_elements, \
                                self.num_steps, self.heat, self.advance_method, sigma)
        outfile = open(filename_referanse, "w")
        # time stepping
        for n in range(0, self.num_steps ):
            self.n = n
            #if n > 6:
            #    break
            if n == 0 and advance_method == "FE":
                # loser for u_e i forste tidssteg med conjugate gradient metoden
                solve(A2, u.vector(), -1*A1*self.v_1_vector, "cg")
                self.u_1_vector = u.vector()
                print self.u_1_vector.array()
            v.vector()[:], w.vector()[:], u.vector()[:] = self.advance()
            self.w_1_vector = w.vector()
            self.v_1_vector = v.vector()
            self.u_1_vector = u.vector()
            self.tid += self.dt
            self.exact_expression.t += self.dt
            if n%20 == 0:
                  print self.tid
            # break if numerical solution diverges
            if abs(np.sum(self.v_1_vector.array())/len(self.v_1_vector.array())) > 10:
                print "break"
                break
            if self.plotting and n%1 == 0:
                plot(v)
        #vtk_v << v
        #write numerical solution to file
        outfile.write("numerical solution at t=%1.8f, sigma = %1.8f" % (self.tid, sigma))
        outfile.write("\n")
        s = "%14.8f" % v.vector().array()[0]
        for j in range(1,len(v.vector().array())):
            s += "%14.8f " % v.vector().array()[j]
        outfile.write(s + "\n")
        outfile.close()


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
        (psi1, psi2) = TestFunctions(V)
        v_var, u_var = split(u_v)
        F = (v_var - _v)/dt*psi1*dx + inner(M_i*grad(v_var),grad(psi1))*dx + inner(M_i*grad(u_var),grad(psi1))*dx + \
            inner(M_i*grad(v_var), grad(psi2))*dx + inner((M_i + M_e)*grad(u_var), grad(psi2))*dx
        #u_var, v_var = TrialFunctions(V)             # her erstattes v- og u-funksjoner som er definert ovenfor
        #elliptic = inner(M_i*grad(v_var),grad(psi1))*dx + inner(M_i*grad(u_var),grad(psi1))*dx
        #parabolic = inner(M_i*grad(v_var), grad(psi2))*dx + inner((M_i + M_e)*grad(u_var), grad(psi2))*dx
        #G = (v_var - _v)*psi1*dx + dt*elliptic + dt*parabolic
        #a, L = system(G)
        #pde = LinearVariationalProblem(a, L, u_v)
        #solver = LinearVariationalSolver(pde)
        #solver.solve()
        #print u_v.vector().array()
        solve(F == 0, u_v)
        u_v_array = u_v.vector().array()
        #_v_array = u_v_array[range(self.n%2 == 0,len(u_v_array),2)]
        #_v.vector()[:] =  _v_array[range(len(_v_array)-1, -1, -1)]
        #_u_array = u_v_array[range(self.n%2 == 0,len(u_v_array),2)]
        #_u.vector()[:] =  _u_array[range(len(_u_array)-1, -1, -1)]
        v_split, u_split = u_v.split(deepcopy=True)
        return v_split.vector(), w_n_s1, u_split.vector()


vtk_v = File("v_bidomain.pvd")

# create Expressions
if two_d:
    if heat:
        I_v_expression = Expression("sigma*cos(pi*x[0])", degree = 2, sigma = sigma)
    else:
        I_v_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = 2)
    I_w_expression = Expression("x[0] < 0.2 && x[1] < 0.2 ? 0: 0", degree = 1)
    I_u_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = 2)
    mesh = UnitSquareMesh(Nx, Ny)

else:
    I_v_expression = Expression("x[0] < 0.2 ? 1: 0, 0", degree = 1)
    I_w_expression = Expression("x[0] < 0.2 ? 0: 0", degree = 2)
    I_u_expression = Expression("x[0] < 0.2 ? 1: 0", degree = 2)
    mesh = UnitIntervalMesh(Nx)

exact_expression = Expression("exp(-pi*pi*t)*cos(pi*x[0])", degree = 2, t = 0)


#create function space
mesh_mix = UnitIntervalMesh(Nx)
V_fe_mix = FiniteElement("CG", mesh.ufl_cell(), degree)
V = FunctionSpace(mesh_mix, MixedElement(V_fe_mix, V_fe_mix))
W = FunctionSpace(mesh_mix, V_fe_mix)
w = Function(W)
u = Function(W)
v = Function(W)
_v = Function(W)
_u = Function(W)
u_v = Function(V)

#interpolate initial condition
v_0 = interpolate(I_v_expression, W)
w_0 = interpolate(I_w_expression, W)
u_0 = interpolate(I_u_expression, W)


many_object_os = OperatorSplitting(T = T, num_steps = num_steps, plotting = plotting, ref = True,\
                                 num_elements = Nx,system_ = system_, heat = heat, exact_expression = exact_expression,
                                 advance_method = advance_method)

if advance_method == "FE":
    many_solution_fe = many_object_fe.solver()
else:
    many_solution_os = many_object_os.solver()

end = time.time()
print end - start
interactive()


"""

# define variational problem for implicit scheme
parabolic = inner(M_i*grad(v), grad(psi1))*dx + inner(M_i*grad(u_e_os), grad(psi1))*dz
elliptic = inner(M_i*grad(v), grad(psi2))*dz + inner((M_i + M_e)*grad(u_e_os), grad(psi2))*dz
G = parabolic + elliptic

block_mat = BlockMatrix(2, 2)
block_mat[0,0] = M - dt*A1
block_mat[0,1] = -dt*A1
block_mat[1,0] = A1
block_mat[1,1] = A2
block_vector = BlockVector(2)
"""
