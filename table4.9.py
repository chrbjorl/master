from conv_rates import *

two_d = False
mono = False
advance_method = "OS"
num_steps = [10*4**i for i in range(4)]
num_elements = [510*2**i for i in range(4)]

table(num_steps = num_steps, num_elements = num_elements, advance_method = advance_method,
    two_d = two_d, mono = mono, plotting = False)
