from fenics import *
import numpy as np
import time
import numpy.linalg as LA
from scipy.linalg import solve_banded
import scipy.sparse.linalg
import scipy.sparse
import dolfin
import ufl
sigma = 1
c_1 = 200
a_1 = 0.1
c_2 = 200
c_3 = 1
b = 1

level = 6
Nx = 20
Ny = 20
two_d = True
degree = 2 if two_d else 1
print degree
T = 1
num_steps = 100000
l2_step = int(num_steps/2.)
system = True
heat = False
plotting = False
advance_method = "FE"
dt = T/float(num_steps)

def scale_m(n):
    return 80 + 100*(n - 2) + 40*((n + 3)*(n - 2)/2 - 3*(n - 2))

class ODEsolver:
    def __init__(self, T, num_steps, plotting, ref, l2_step, num_elements,
                system, heat, exact_expression, advance_method):
        self.T, self.num_steps = T, num_steps
        self.plotting, self.ref = plotting, ref
        self.v_1_vector = v_0.vector()
        self.w_1_vector = w_0.vector()
        self.dt = T/float(num_steps)
        self.l2_step = l2_step
        self.l2_time = l2_step/float(num_steps)*T
        self.num_elements = num_elements
        self.tid = 0
        self.system = system
        self.heat = heat
        self.exact_expression = exact_expression
        self.advance_method = advance_method


    def solver(self):
        ref = self.ref
        filename_referanse = "%sd_%s_%s_%1.9f_%s_%s.dat" % (degree, self.num_elements, \
                                self.num_steps, self.l2_time, self.heat, self.advance_method)
        outfile = open(filename_referanse, "w")
        for n in range(0, self.l2_step + 2):
            if system:
                v.vector()[:], w.vector()[:] = self.advance()
                self.w_1_vector = w.vector()
            else:
                v.vector()[:] = self.advance()
            self.v_1_vector = v.vector()
            self.tid += self.dt
            self.exact_expression.t += self.dt
            if n%10 == 0:
                print self.tid
            if abs(np.sum(self.v_1_vector.array())/len(self.v_1_vector.array())) > 10000:
                print "hei"
                break
            if n == self.l2_step - 1:
                u_e = project(self.exact_expression, V)
                #write numerical solution to file
                outfile.write("numerical solution at t=%1.8f" % self.tid)
                outfile.write("\n")
                s = "%14.8f" % v.vector().array()[0]
                for j in range(1,len(v.vector().array())):
                    s += "%14.8f" % v.vector().array()[j]
                outfile.write(s + "\n")
                if self.heat:
                    # write exact solution to file
                    outfile.write("exact solution at t=%1.8f" % self.tid)
                    outfile.write("\n")
                    s = "%14.8f" % u_e.vector().array()[0]
                    for j in range(1,len(u_e.vector().array())):
                        s += "%14.8f" % u_e.vector().array()[j]
                    outfile.write(s + "\n")
                    outfile.close()
                break
            if self.plotting and n%10 == 0:
                plot(v)

class OperatorSplitting(ODEsolver):
    """ Bruker Godunov-splitting og baklengs Euler
    """

    def advance(self):
        v_1_vector = self.v_1_vector
        w_1_vector = self.w_1_vector
        dt = self.dt
        # step 1
        RHS = M*v_1_vector
        solve(LHS, v_n_s2.vector(), RHS)
        if self.system:
            # step 2
            v_n_s3 = v_n_s2.vector() + dt*(c_1*v_n_s2.vector()*\
                (v_n_s2.vector() - a_1)*(1 - v_n_s2.vector()) - c_2*v_n_s2.vector()*w_1_vector)

            w_n_s3 = dt*b*(v_n_s2.vector() - c_3*w_1_vector) + w_1_vector
            return v_n_s3, w_n_s3
        else:
            return v_n_s2.vector()

class ForwardEuler(ODEsolver):
    def advance(self):
        if self.system:
            w_1_vector = self.w_1_vector
        v_1_vector = self.v_1_vector
        dt = self.dt
        if self.system:
            product1 = dt*A*v_1_vector
            return M_inv*product1 + v_1_vector + dt*c_1*v_1_vector*\
            (v_1_vector - a_1)*(1 - v_1_vector) - dt*c_2*v_1_vector*w_1_vector, \
            dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector
        if self.heat:
            product1 = dt*A*v_1_vector
            return M_inv*product1 + v_1_vector
        if not self.system and not self.heat:
            return scale*20*dt*A*v_1_vector + v_1_vector + dt*c_1*v_1_vector*\
            (v_1_vector - a_1)*(1 - v_1_vector)

def l2(x):
    return np.sqrt(np.sum(x**2)/float(len(x)))

# create Expressions
#I_v_expression = Expression("cos(pi*x[0])", degree = 2)
if two_d:
    if heat:
        I_v_expression = Expression("cos(pi*x[0])", degree = 2)
    else:
        I_v_expression = Expression("x[0] < 0.5 && x[1] < 0.5 ? 1: 0", degree = 2)
    I_w_expression = Expression("x[0] < 0.5 && x[1] < 0.5 ? 0: 0", degree = 1)
    mesh = UnitSquareMesh(Nx, Ny)

else:
    if heat:
        I_v_expression = Expression("cos(pi*x[0])", degree = 2)
    else:
        I_v_expression = Expression("x[0] < 0.5 ? 1: 0", degree = 2)
    I_w_expression = Expression("x[0] < 0.5 ? 0: 0", degree = 2)
    mesh = UnitIntervalMesh(Nx)

exact_expression = Expression("exp(-pi*pi*t)*cos(pi*x[0])", degree = 2, t = 0)

#create function space
V = FunctionSpace(mesh, "P", degree)

#interpolate initial condition
v_0 = interpolate(I_v_expression, V)
w_0 = interpolate(I_w_expression, V)

# Define variational problem
v = TrialFunction(V)
psi = TestFunction(V)
v_n_s2 = Function(V)
a = dot(-sigma*grad(v), grad(psi))*dx
m = dot(v, psi)*dx
# assemble A outside time loop, since A is time-independent
A = assemble(a)
M = assemble(m)
v = Function(V)
w = Function(V)

q = TrialFunction(V)
z = TestFunction(V)
form = inner(q,z)*dx
I = assemble(form)
I.set_diagonal(interpolate(Constant(1), V).vector())
LHS = M - A*dt

diag = [M.array()[j][j] for j in range(len(M.array()))]
diag2 = [sum(M.array()[j]) for j in range(len(M.array()))]
c = sum(diag2)/sum(diag)
diag = np.asarray(diag)
diag *= c
diag = diag**(-1)
v.vector()[:] = diag[:]
M_inv = M
M_inv.zero()
M_inv.set_diagonal(v.vector())
M = assemble(m)

start = time.time()

many_object_fe = ForwardEuler(T = T, num_steps = num_steps, plotting = plotting, ref = True,\
                            l2_step = l2_step, num_elements = Nx,
                           system = system, heat = heat, exact_expression = exact_expression,
                            advance_method = advance_method)

many_object_os = OperatorSplitting(T = T, num_steps = num_steps, plotting = plotting, ref = True,\
                                l2_step = l2_step, num_elements = Nx,
                                system = system, heat = heat, exact_expression = exact_expression,
                                advance_method = advance_method)

if advance_method == "FE":
    many_solution_fe = many_object_fe.solver()
else:
    many_solution_os = many_object_os.solver()

end = time.time()
print end - start

#interactive()
