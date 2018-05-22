from conv_rates import *

two_d = 0
mono = False
advance_method = "FE"
jacobi_ = True
num_jacobi_iter = [1, 1, 1, 1]
num_steps = [425*2**i for i in range(4)]
num_elements = [50*2**i for i in range(4)]

table(num_steps = num_steps, num_elements = num_elements, advance_method = advance_method,
    two_d = False, mono = mono, plotting = False, jacobi_ = jacobi_,
    num_jacobi_iter = num_jacobi_iter)
