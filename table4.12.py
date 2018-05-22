from conv_rates import *

two_d = 1
mono = False
advance_method = "OS"
num_steps = [8*4**i for i in range(3)]
num_elements = [60*2**i for i in range(3)]

table(num_steps = num_steps, num_elements = num_elements, advance_method = advance_method,
    two_d = two_d, mono = mono, plotting = False)
