from fenics import *
import numpy as np
import time
import matplotlib.pyplot as plt
import os

M_i  = 3
M_e  = 3
c_1 = 200
a_1 = 0.1
c_2 = 200
c_3 = 1
b = 1
num_jacobi_iter = 1

Nx = 20
Ny = Nx

degree = 1
two_d = 1
T = 0.075
num_steps = 10000

plotting = True
advance_method = "FE"           #"FE" if ForwardEuler and "OS" if OperatorSplitting
dt = T/float(num_steps)
vtk_v = File("v_bidomain.pvd")

def jacobi(u_0, K, v_fe):
    u_old = u_0
    product1 = A1*v_fe
    for k in range(K):
        u_new = -1*D_inv*(R*u_old - product1)
        # u_new = u_old - D_inv*(A2*u_old - product1)
        u_old = u_new
    return u_new

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
            if n%200 == 0:
                  print self.tid
            # break if numerical solution diverges
            if abs(np.sum(self.v_1_vector.get_local())/len(self.v_1_vector.get_local())) > 10:
                print "break"
                break
            if self.plotting and n%int(num_steps/4) == 0:
                plot(v)
                plt.show()
        vtk_v << v
        #write numerical solution to file
        output_file = XDMFFile(filename + ".xdmf")
        output_file.write_checkpoint(v,'v')
        File(filename + '_mesh.xml.gz') << v.function_space().mesh()


class ForwardEuler(ODEsolver):
    def advance(self):
        w_1_vector = self.w_1_vector
        v_1_vector = self.v_1_vector
        u_1_vector = self.u_1_vector
        dt = self.dt
        product1 = dt*A1*(v_1_vector + u_1_vector)
        v_fe = M_inv*product1 + v_1_vector + dt*c_1*v_1_vector*\
                (v_1_vector - a_1)*(1 - v_1_vector) - dt*c_2*v_1_vector*w_1_vector
        w_fe = dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector
        # use one Jacobi iteration with the solution at the previous time step as
        # initial guess
        u_fe = jacobi(u_1_vector, num_jacobi_iter, v_fe)
        return v_fe, w_fe, u_fe


# create Expressionsf
if two_d:
    I_v_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = 2)
    I_w_expression = Expression("x[0] < 0.2 && x[1] < 0.2 ? 0: 0", degree = 1)
    I_u_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = 2)
    mesh = UnitSquareMesh(Nx, Ny)

else:
    I_v_expression = Expression("x[0] < 0.2 ? 1: 0", degree = 2)
    I_w_expression = Expression("x[0] < 0.2 ? 0: 0", degree = 2)
    I_u_expression = Expression("x[0] < 0.2 ? 1: 0", degree = 2)
    mesh = UnitIntervalMesh(Nx)
start = time.time()
V = FunctionSpace(mesh, "P", degree)

#interpolate initial condition
v_0 = interpolate(I_v_expression, V)
w_0 = interpolate(I_w_expression, V)
u_0 = interpolate(I_u_expression, V)

# Define variational problem used for explicit scheme
v = TrialFunction(V)
u_e = TrialFunction(V)
psi = TestFunction(V)

a1 = dot(-M_i*grad(v), grad(psi))*dx
m = dot(v, psi)*dx
a2 = dot((M_i + M_e)*grad(u_e), grad(psi))*dx

# assemble A outside time loop, since A is time-independent
A1 = assemble(a1)
A2 = assemble(a2)
M = assemble(m)
v = Function(V)
w = Function(V)
y = Function(V)
u = Function(V)
zeros = Function(V)
LHS = M - A1*dt

# lumped mass matrix
# diagonal elements of M
M.get_diagonal(y.vector())
diag =  y.vector().array()
#create identity matrix
I = M
I.zero()
I.set_diagonal(interpolate(Constant(1), V).vector())
I.get_diagonal(y.vector())
# diag2 contains the row sums of M
M = assemble(m)
S = M*y.vector()
diag2 =  S.array()
c = sum(diag2)/sum(diag)
diag *= c
diag = diag**(-1)
v.vector()[:] = diag[:]
M_inv = I
M_inv.zero()
M_inv.set_diagonal(v.vector())

#D_inv er en diagonalmatrise som brukes i  Jacobi-iterasjoner i det eksplisitte skjemaet
A2.get_diagonal(y.vector())
diag3 =  y.vector().array()
diag3 = diag3**(-1)
D_inv = A2.copy()
D_inv.zero()
v.vector()[:] = diag3[:]
D_inv.set_diagonal(v.vector())
D_inv.get_diagonal(y.vector())

# R er en matrise som brukes i Jacobi-iterasjonene
R = A2.copy()
zeros.vector()[:] = np.zeros(len(zeros.vector().array()))
R.set_diagonal(zeros.vector())
object_fe = ForwardEuler(T = T, num_steps = num_steps, plotting = plotting,
                             num_elements = Nx, advance_method = advance_method)

solution_fe = object_fe.solver()

end = time.time()
print end - start
#interactive()
