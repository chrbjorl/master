from fenics import *
import numpy as np
import time
import matplotlib.pyplot as plt

sigma = 0.03
#c_1 = 200
#a_1 = 0.1
#c_2 = 200
#c_3 = 1
#b = 1

Nx = 40
Ny = Nx

degree = 1
two_d = 0
T = 0.2
num_steps = 1000

system = False         # fitzhugh-nagumo if system = True
heat = False
plotting = True
advance_method = "FE"           #"FE" if ForwardEuler and "OS" if OperatorSplitting
dt = T/float(num_steps)

class ODEsolver:
    def __init__(self, T, num_steps, plotting, num_elements,
                system, heat, advance_method, two_d, sigma,
                c_1 = 200, a_1 = 0.1, c_2 = 200, c_3 = 1, b = 1):
        self.T, self.num_steps = T, num_steps
        self.plotting = plotting
        self.dt = T/float(num_steps)
        self.num_elements = num_elements
        self.tid = 0
        self.system = system
        self.heat = heat
        self.advance_method = advance_method
        self.two_d = two_d
        self.c_1, self.a_1, self.c_2, self.c_3, self.b = c_1, a_1, c_2, c_3, b
        Nx = num_elements
        Ny = Nx

        # create Expressions
        if two_d:
            if heat:
                I_v_expression = Expression("sigma*cos(pi*x[0])", degree = 2, sigma = sigma)
            else:
                I_v_expression = Expression("pow(x[0], 2) + pow(x[1], 2) < 0.2 ? 1: 0", degree = 2)
            I_w_expression = Expression("x[0] < 0.2 && x[1] < 0.2 ? 0: 0", degree = 1)
            mesh = UnitSquareMesh(Nx, Ny)

        else:
            if heat:
                I_v_expression = Expression("sigma*cos(pi*x[0])", degree = 2, sigma = sigma)
            else:
                I_v_expression = Expression("x[0] < 0.2 ? 1: 0", degree = 2)
            I_w_expression = Expression("x[0] < 0.2 ? 0: 0", degree = 2)
            mesh = UnitIntervalMesh(Nx)

        exact_expression = Expression("sigma*exp(-pi*pi*t)*cos(pi*x[0])", degree = 2, t = 0,
                                        sigma = sigma)

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
        y = Function(V)
        z = Function(V)
        LHS = M - A*dt

        # lumped mass matrix
        # diagonal elements of M
        M.get_diagonal(y.vector())
        diag =  y.vector().get_local()
        #create identity matrix
        I = M.copy()
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
        # henter ut diagonalen til lumped mass matrix
        M_inv.get_diagonal(z.vector())
        M_inv_diag_vector = z.vector()
        self.v_1_vector = v_0.vector()
        self.w_1_vector = w_0.vector()
        self.M_inv = M_inv
        self.M = M
        self.A = A
        self.LHS = LHS
        self.sigma = sigma
        self.v, self.w, self.V = v, w, V
        self.exact_expression = exact_expression



    def solver(self):
        v, w = self.v, self.w
        filename = "%sd_%s_%s_%s_%s_%s_%1.3f" % (2**self.two_d, self.num_elements, \
                                self.num_steps, self.heat, self.system, self.advance_method, self.sigma)
        outfile = open(filename, "w")
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
            if n%100 == 0:
                  print self.tid
            # break if numerical solution diverges
            if abs(np.sum(self.v_1_vector.get_local())/len(self.v_1_vector.get_local())) > 2:
                print "break"
                break
            if self.plotting and n%int(self.num_steps/2) == 0:
                plot(v)
                plt.show()
        u_e = project(self.exact_expression, self.V)
        #write numerical solution to file
        output_file = XDMFFile(filename + ".xdmf")
        output_file.write_checkpoint(v,'v')
        File(filename + '_mesh.xml.gz') << v.function_space().mesh()
        if self.heat:
            # write exact solution to file
            outfile.write("exact solution at t=%1.8f" % self.tid)
            outfile.write("\n")
            s = "%14.8f" % u_e.vector().get_local()[0]
            for j in range(1,len(u_e.vector().get_local())):
                s += "%14.8f" % u_e.vector().get_local()[j]
            outfile.write(s + "\n")
            outfile.close()

class OperatorSplitting(ODEsolver):
    """ Bruker Godunov-splitting og baklengs Euler
    """

    def advance(self):
        v_1_vector = self.v_1_vector
        w_1_vector = self.w_1_vector
        dt = self.dt
        LHS = self.LHS
        M = self.M
        # step 1
        RHS = M*v_1_vector
        solve(LHS, v_n_s2.vector(), RHS, "cg")
        c_1, a_1, c_2, c_3, b = self.c_1, self.a_1, self.c_2, self.c_3, self.b
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
        dt = self.dt
        LHS = self.LHS
        M = self.M
        M_inv = self.M_inv
        A = self.A
        product1 = dt*A*v_1_vector
        c_1, a_1, c_2, c_3, b = self.c_1, self.a_1, self.c_2, self.c_3, self.b
        if self.system:
            return M_inv*product1 + v_1_vector + dt*c_1*v_1_vector*\
            (v_1_vector - a_1)*(1 - v_1_vector) - dt*c_2*v_1_vector*w_1_vector, \
            dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector
        if self.heat:
            product1 = dt*A*v_1_vector
            return M_inv*product1 + v_1_vector
        if not self.heat and not self.system:
            return M_inv*product1 + v_1_vector + dt*c_1*v_1_vector*\
                    (v_1_vector - a_1)*(1 - v_1_vector)


class ForwardEulerSystem(ODEsolver):
    def advance(self):
        w_1_vector = self.w_1_vector
        v_1_vector = self.v_1_vector
        dt = self.dt
        product1 = dt*A*v_1_vector
        return M_inv*product1 + v_1_vector + dt*c_1*v_1_vector*\
                (v_1_vector - a_1)*(1 - v_1_vector) - dt*c_2*v_1_vector*w_1_vector, \
                dt*b*(v_1_vector - c_3*w_1_vector) + w_1_vector


if __name__ == "__main__":
    start = time.time()
    object_fe = ForwardEuler(T = T, num_steps = num_steps, plotting = plotting, \
                            num_elements = Nx, system = system, heat = heat,
                              advance_method = advance_method, two_d = two_d,
                              sigma = sigma)

    object_os = OperatorSplitting(T = T, num_steps = num_steps, plotting = plotting,
                                    num_elements = Nx, system = system, heat = heat,
                                     advance_method = advance_method, two_d = two_d,
                                     sigma = sigma)
    object_fe_system = ForwardEulerSystem(T = T, num_steps = num_steps, plotting = plotting,
                                    num_elements = Nx, system = system, heat = heat,
                                    advance_method = advance_method, two_d = two_d,
                                    sigma = sigma)
    if advance_method == "FE":
        solution_fe = object_fe.solver()
    elif advance_method == "OS":
        solution_os = object_os.solver()
    elif advance_method == "FE_system":
        solution_os = object_fe_system.solver()

    end = time.time()
    print end - start
    #interactive()
