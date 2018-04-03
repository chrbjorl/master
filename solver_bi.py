from fenics import *
import numpy as np
import time
import scipy.sparse

M_i  = 0.8
M_e  = 0.8
sigma = 1
c_1 = 200
a_1 = 0.1
c_2 = 200
c_3 = 1
b = 1

Nx = 30
Ny = Nx

degree = 1
two_d = 0
T = 0.5
num_steps = 2000

system = True            # fitzhugh-nagumo if system = True
heat = False
plotting = True
advance_method = "FE"           #"FE" if ForwardEuler and "OS" if OperatorSplitting
dt = T/float(num_steps)

class ODEsolver:
    def __init__(self, T, num_steps, plotting, ref, num_elements,
                system, heat, exact_expression, advance_method):
        self.T, self.num_steps = T, num_steps
        self.plotting, self.ref = plotting, ref
        self.v_1_vector = v_0.vector()
        self.w_1_vector = w_0.vector()
        self.dt = T/float(num_steps)
        self.num_elements = num_elements
        self.tid = 0
        self.system = system
        self.heat = heat
        self.exact_expression = exact_expression
        self.advance_method = advance_method


    def solver(self):
        ref = self.ref
        filename_referanse = "%sd_%s_%s_%s_%s_%1.3f.dat" % (2**two_d, self.num_elements, \
                                self.num_steps, self.heat, self.advance_method, sigma)
        outfile = open(filename_referanse, "w")
        # time stepping
        for n in range(0, self.num_steps ):
            if system:
                v.vector()[:], w.vector()[:] = self.advance()
                self.w_1_vector = w.vector()
            else:
                v.vector()[:] = self.advance()
            self.v_1_vector = v.vector()
            self.tid += self.dt
            self.exact_expression.t += self.dt
            if n%20 == 0:
                  print self.tid
            # break if numerical solution diverges
            if abs(np.sum(self.v_1_vector.array())/len(self.v_1_vector.array())) > 2:
                print "break"
                break
            if self.plotting and n%2 == 0:
                plot(v)
        u_e = project(self.exact_expression, V)
        #write numerical solution to file
        outfile.write("numerical solution at t=%1.8f, sigma = %1.8f" % (self.tid, sigma))
        outfile.write("\n")
        s = "%14.8f" % v.vector().array()[0]
        for j in range(1,len(v.vector().array())):
            s += "%14.8f " % v.vector().array()[j]
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

class OperatorSplitting(ODEsolver):
    """ Bruker Godunov-splitting og baklengs Euler
    """

    def advance(self):
        v_1_vector = self.v_1_vector
        w_1_vector = self.w_1_vector
        dt = self.dt
        # step 1
        RHS = M*v_1_vector
        solve(LHS, v_n_s2.vector(), RHS, "cg")
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
        w_1_vector = self.w_1_vector
        v_1_vector = self.v_1_vector
        u_1_vector = self.u_1_vector
        dt = self.dt
        product1 = dt*A1*(v_1_vector + u_1_vector)
        v_fe = M_inv_diag_vector*product1 + v_1_vector + dt*c_1*v_1_vector*\
                (v_1_vector - a_1)*(1 - v_1_vector) - dt*c_2*v_1_vector*w_1_vector
        w_fe = dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector
        # use one Jacobi iteration with the solution at the previous time step as
        # initial guess
        u_fe = u_1_vector - D1*(A2*u_1_vector + A1*v_fe)
        return v_fe, w_fe, u_fe

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
    if heat:
        I_v_expression = Expression("sigma*cos(pi*x[0])", degree = 2, sigma = sigma)
    else:
        I_v_expression = Expression("x[0] < 0.2 ? 1: 0", degree = 2)
    I_w_expression = Expression("x[0] < 0.2 ? 0: 0", degree = 2)
    I_u_expression = Expression("x[0] < 0.2 ? 1: 0", degree = 2)
    mesh = UnitIntervalMesh(Nx)

exact_expression = Expression("exp(-pi*pi*t)*cos(pi*x[0])", degree = 2, t = 0)

#create function space
V = FunctionSpace(mesh, "P", degree)

#interpolate initial condition
v_0 = interpolate(I_v_expression, V)
w_0 = interpolate(I_w_expression, V)
u_0 = interpolate(I_u_expression, V)

# Define variational problem
v = TrialFunction(V)
psi = TestFunction(V)
v_n_s2 = Function(V)
a1 = dot(-M_i*grad(v), grad(psi))*dx
m = dot(v, psi)*dx
a2 = dot((M_i + M_e)*grad(v), grad(psi))*dx

# assemble A outside time loop, since A is time-independent
A1 = assemble(a1)
A2 = assemble(a2)
M = assemble(m)
v = Function(V)
w = Function(V)
y = Function(V)
u = Function(V)
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

# henter ut diagonalen til massematrisen
M_inv.get_diagonal(z.vector())
print z.vector().array()
M_inv_diag_vector = z.vector()

#D1 er en diagonalmatrise som brukes i  Jacobi-iterasjoner i det eksplisitte skjemaet
A2.get_diagonal(y.vector())
diag3 =  y.vector().array()
diag3 = diag3**(-1)
D1 = A2.copy()
D1.zero()
v.vector()[:] = diag3[:]
D1.set_diagonal(v.vector())
start = time.time()

# sjekker om Jacobi-iterasjonene i det eksplisitte skjemaet vil konvergere
D1_mat = as_backend_type(D1).mat()
A2_mat = as_backend_type(A2).mat()
C = D1_mat.matMult(A2_mat)
C = Matrix(PETScMatrix(C))
print C.array()

many_object_fe = ForwardEuler(T = T, num_steps = num_steps, plotting = plotting, ref = True,\
                             num_elements = Nx, system = system, heat = heat, exact_expression = exact_expression,
                             advance_method = advance_method)

many_object_os = OperatorSplitting(T = T, num_steps = num_steps, plotting = plotting, ref = True,\
                                 num_elements = Nx,system = system, heat = heat, exact_expression = exact_expression,
                                 advance_method = advance_method)

if advance_method == "FE":
    many_solution_fe = many_object_fe.solver()
else:
    many_solution_os = many_object_os.solver()

end = time.time()
print end - start
interactive()
