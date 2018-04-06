from fenics import *
import numpy as np
import time
import scipy.sparse
import block_mat

M_i  = 3
M_e  = 3
sigma = 1
c_1 = 200
a_1 = 0.1
c_2 = 200
c_3 = 1
b = 1
num_jacobi_iter = 2

Nx = 10
Ny = Nx

degree = 1
two_d = 0
T = 0.2
num_steps = 1000

system = True            # fitzhugh-nagumo if system = True
heat = False
plotting = True
advance_method = "OS"           #"FE" if ForwardEuler and "OS" if OperatorSplitting
dt = T/float(num_steps)

def jacobi(u_0, K, v_fe):
    u_old = u_0
    product1 = A1*v_fe
    for k in range(K):
        #product1 = -1*A1*v_fe
        #u_new = -0.01041667*(product1 + R*u_old)
        u_new = u_old + alfa_jacobi*(product1 - A2*u_old)
        u_old = u_new
    return u_new

class ODEsolver:
    def __init__(self, T, num_steps, plotting, ref, num_elements,
                system, heat, exact_expression, advance_method):
        self.T, self.num_steps = T, num_steps
        self.plotting, self.ref = plotting, ref
        self.v_1_vector = v_0.vector()
        self.w_1_vector = w_0.vector()
        self.u_1_vector = u_0.vector()
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
            self.n = n
            #if n > 6:
            #    break
            if n == 0:
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
            if self.plotting and n%10 == 0:
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
        v_n_s1 = v_1_vector + dt*(c_1*v_1_vector*\
                (v_1_vector - a_1)*(1 - v_1_vector) - c_2*v_1_vector*w_1_vector)

        w_n_s1 = dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector
        # step 2
        block_vector[0] = M*v_n_s1
        block_vector[1] = zeros.vector()
        RHS = M*v_1_vector
        print len(x.vector().array())
        print len(v_n_s1.array())
        solve(LHS, x.vector(), RHS, "cg")
        return v_n_s3, w_n_s3, x.vector()

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
        solve(A2, y.vector(), -1*A1*v_fe, "cg")
        #vec1 = A2*u_fe
        #vec2 = -1*A1*v_fe
        return v_fe, w_fe, y.vector()


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

mesh_x = UnitIntervalMesh(2*Nx + 1)
V_x = FunctionSpace(mesh_x, "P", degree)
x = Function(V_x)
#interpolate initial condition
v_0 = interpolate(I_v_expression, V)
w_0 = interpolate(I_w_expression, V)
u_0 = interpolate(I_u_expression, V)

# Define variational problem
v = TrialFunction(V)
psi = TestFunction(V)
u_e_os = Function(V)
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


#D1 er en diagonalmatrise som brukes i  Jacobi-iterasjoner i det eksplisitte skjemaet
A2.get_diagonal(y.vector())
diag3 =  y.vector().array()
diag3 = diag3**(-1)
D1 = A2.copy()
D1.zero()
v.vector()[:] = diag3[:]
D1.set_diagonal(v.vector())
D1.get_diagonal(y.vector())
alfa_jacobi =  y.vector().array()
alfa_jacobi = alfa_jacobi[5]

print "alfa jacobi: %s" % alfa_jacobi

# R er en matrise som brukes i Jacobi-iterasjonene
R = A2.copy()
zeros.vector()[:] = np.zeros(len(zeros.vector().array()))
R.set_diagonal(zeros.vector())
start = time.time()

# sjekker om Jacobi-iterasjonene i det eksplisitte skjemaet vil konvergere
D1_mat = as_backend_type(D1).mat()
A2_mat = as_backend_type(A2).mat()
R_mat = as_backend_type(R).mat()
C1 = D1_mat.matMult(A2_mat)
C1 = Matrix(PETScMatrix(C1))
C2 = PETScMatrix(D1_mat.matMult(R_mat))


#C2 = Matrix(PETScMatrix(C2))


eigensolver = SLEPcEigenSolver(C2)
eigensolver.solve()
r, c, rx, cx = eigensolver.get_eigenpair(1)
print "largest eigenvalue of g: %1.7f" % r
block_mat = BlockMatrix(2, 2)
block_mat[0,0] = M - dt*A1
block_mat[0,1] = -dt*A1
block_mat[1,0] = A1
block_mat[1,1] = A2
block_vector = BlockVector(2)


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
