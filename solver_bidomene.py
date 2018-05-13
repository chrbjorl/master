from fenics import *
import numpy as np
import time
import matplotlib.pyplot as plt
import os


class PDEsolver:
    def __init__(self, T, num_steps, plotting_v, num_elements,
                advance_method, two_d, degree, num_jacobi_iter, jacobi_, M_i = 0.1, M_e = 0.1,
                c_1 = 200, a_1 = 0.1, c_2 = 200, c_3 = 1, b = 1):
        self.T, self.num_steps = T, num_steps
        self.plotting_v = plotting_v
        self.dt = T/float(num_steps)
        self.num_elements = num_elements
        self.tid = 0
        self.advance_method = advance_method
        self.two_d = two_d
        self.c_1, self.a_1, self.c_2, self.c_3, self.b = c_1, a_1, c_2, c_3, b
        self.start = time.time()
        Nx = num_elements
        Ny = Nx
        self.start = time.time()

        # create Expressions
        if two_d:
            I_v_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = degree)
            I_w_expression = Expression("x[0] < 0.2 && x[1] < 0.2 ? 0: 0", degree = degree)
            I_u_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = degree)
            mesh = UnitSquareMesh(Nx, Ny)

        else:
            I_v_expression = Expression("x[0] < 0.2 ? 1: 0", degree = degree)
            I_w_expression = Expression("x[0] < 0.2 ? 0: 0", degree = degree)
            I_u_expression = Expression("x[0] < 0.2 ? 1: 0", degree = degree)
            mesh = UnitIntervalMesh(Nx)

        #interpolate initial condition

        if self.advance_method == "FE":# Define variational problem used for explicit scheme
            V = FunctionSpace(mesh, "P", degree)
            v = TrialFunction(V)
            u_e = TrialFunction(V)
            psi = TestFunction(V)

            v_0 = interpolate(I_v_expression, V)
            w_0 = interpolate(I_w_expression, V)
            u_0 = interpolate(I_u_expression, V)

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

            # lumped mass matrix
            # diagonal elements of M
            M.get_diagonal(y.vector())
            diag =  y.vector().get_local()
            #create identity matrix
            I = M
            I.zero()
            I.set_diagonal(interpolate(Constant(1), V).vector())
            I.get_diagonal(y.vector())
            # diag2 contains the row sums of M
            M = assemble(m)
            S = M*y.vector()
            diag2 =  S.get_local()
            c = sum(diag2)/sum(diag)
            diag *= c
            diag = diag**(-1)
            v.vector()[:] = diag[:]
            M_inv = I
            M_inv.zero()
            M_inv.set_diagonal(v.vector())

            #D_inv er en diagonalmatrise som brukes i  Jacobi-iterasjoner i det eksplisitte skjemaet
            A2.get_diagonal(y.vector())
            diag3 =  y.vector().get_local()
            diag3 = diag3**(-1)
            D_inv = A2.copy()
            D_inv.zero()
            v.vector()[:] = diag3[:]
            D_inv.set_diagonal(v.vector())
            D_inv.get_diagonal(y.vector())

            # R er en matrise som brukes i Jacobi-iterasjonene
            R = A2.copy()
            zeros.vector()[:] = np.zeros(len(zeros.vector().get_local()))
            R.set_diagonal(zeros.vector())
            self.u_1_vector = u_0.vector()

            self.R, self.A2, self.A1, self.D_inv, self.M_inv = R, A2, A1, D_inv, M_inv
            self.num_jacobi_iter, self.jacobi_ = num_jacobi_iter, jacobi_
        else:
            #create function space
            V_fe_mix = FiniteElement("CG", mesh.ufl_cell(), degree)
            V = FunctionSpace(mesh, MixedElement(V_fe_mix, V_fe_mix))
            W = FunctionSpace(mesh, V_fe_mix)
            w = Function(W)
            u = Function(W)
            v = Function(W)
            _v = Function(W)
            _u = Function(W)
            u_v = Function(V)

            #interpolate initial condition
            v_0 = interpolate(I_v_expression, W)
            w_0 = interpolate(I_w_expression, W)
            self._v, self._u, self.u_v, self.W = _v, _u, u_v, W
            self.M_i, self.M_e = M_i, M_e
            u_0 = interpolate(I_u_expression, W)
        self.v_1_vector = v_0.vector()
        self.w_1_vector = w_0.vector()
        self.v, self.w, self.u, self.V = v, w, u, V
        self.num_jacobi_iter, self.jacobi_ = num_jacobi_iter, jacobi_

    def solver(self):
        v, w, u = self.v, self.w, self.u
        if self.advance_method == "FE":
            filename = "bi_%sd_%s_%s_%s_%s_%s" % (2**self.two_d, self.num_elements,
                                self.num_steps, self.advance_method, self.num_jacobi_iter,
                                self.jacobi_)
        else:
            filename = "bi_%sd_%s_%s_%s" % (2**self.two_d, self.num_elements,
                                            self.num_steps, self.advance_method)
        # time stepping
        for n in range(0, self.num_steps ):
            self.n = n
            if n == 0 and self.advance_method == "FE":
                # loser for u_e i forste tidssteg med conjugate gradient metoden
                solve(self.A2, u.vector(), -1*self.A1*self.v_1_vector, "cg")
                self.u_1_vector = u.vector()
                print self.u_1_vector.get_local()
            v.vector()[:], w.vector()[:], u.vector()[:] = self.advance()
            self.w_1_vector = w.vector()
            self.v_1_vector = v.vector()
            self.u_1_vector = u.vector()
            self.tid += self.dt
            if n%10 == 0:
                  print self.tid
            # break if numerical solution diverges
            if abs(np.sum(self.v_1_vector.get_local())/len(self.v_1_vector.get_local())) > 10:
                print "break"
                break
            if self.plotting_v and n%int(num_steps/3) == 0:
                plot(v)
                plt.show()
        end = time.time()
        print end - self.start
        output_file = XDMFFile(filename + ".xdmf")
        output_file.write_checkpoint(v,'v')
        File(filename + '_mesh_v.xml.gz') << v.function_space().mesh()

class ForwardEuler(PDEsolver):
    def advance(self):
        w_1_vector = self.w_1_vector
        v_1_vector = self.v_1_vector
        u_1_vector = self.u_1_vector
        dt = self.dt
        R, A2, A1, D_inv, M_inv = self.R, self.A2, self.A1, self.D_inv, self.M_inv
        c_1, a_1, c_2, c_3, b = self.c_1, self.a_1, self.c_2, self.c_3, self.b
        jacobi_ = self.jacobi_
        product1 = dt*A1*(v_1_vector + u_1_vector)
        v_fe = M_inv*product1 + v_1_vector + dt*c_1*v_1_vector*\
                (v_1_vector - a_1)*(1 - v_1_vector) - dt*c_2*v_1_vector*w_1_vector
        w_fe = dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector
        # use one Jacobi iteration with the solution at the previous time step as
        # initial guess
        if jacobi_:
            product1 = A1*v_fe
            for k in range(self.num_jacobi_iter):
                u_fe = -1*D_inv*(R*u_1_vector - product1)
                #u_new = u_old - D_inv*(A2*u_old - product1)
                u_1_vector = u_fe
                #u_fe = jacobi(u_1_vector, num_jacobi_iter, v_fe)
        else:
            RHS = A1*v_fe
            solve(A2, u.vector(), RHS, "cg")
            u_fe = u.vector()
        return v_fe, w_fe, u_fe

class OperatorSplitting(PDEsolver):
    """ Bruker Godunov-splitting og baklengs Euler
    """
    def advance(self):
        _v, _u, u_v = self._v, self._u, self.u_v
        c_1, a_1, c_2, c_3, b = self.c_1, self.a_1, self.c_2, self.c_3, self.b
        V, W = self.V, self.W
        v_1_vector = self.v_1_vector
        w_1_vector = self.w_1_vector
        dt = self.dt
        M_i, M_e = self.M_i, self.M_e
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
        solve(F == 0, u_v)
        v_split, u_split = u_v.split(deepcopy=True)
        return v_split.vector(), w_n_s1, u_split.vector()

if __name__ == "__main__":
    M_i  = 0.1
    M_e  = 0.1
    c_1 = 200
    a_1 = 0.1
    c_2 = 200
    c_3 = 1
    b = 1
    num_jacobi_iter = 1

    i = 0
    Nx = 50*2**i     #20*2**i
    Ny = Nx

    degree = 1   # ikke endre denne!
    two_d = 1
    T = 0.1
    num_steps = 15*2**i

    jacobi_ = True
    plotting_v = False
    plotting_u = False
    advance_method = "FE"           #"FE" if ForwardEuler and "OS" if OperatorSplitting
    dt = T/float(num_steps)
    object_fe = ForwardEuler(T = T, num_steps = num_steps, plotting_u = plotting_u,
                             plotting_v = plotting_v, num_elements = Nx,
                              advance_method = advance_method)
    object_fe.solver()
