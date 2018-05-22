from conv_rates import *

two_d = 1
mono = False
advance_method = "FE"
jacobi_ = True
num_jacobi_iter = [1,1,1]
num_steps = [65*2**i for i in range(3)]
num_elements = [20*2**i for i in range(3)]

table(num_steps = num_steps, num_elements = num_elements, advance_method = advance_method,
    two_d = two_d, mono = mono, plotting = False, T = 0.1, jacobi_ = jacobi_,
        num_jacobi_iter = num_jacobi_iter)
