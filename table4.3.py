from conv_rates import *

two_d = False
mono = True
advance_method = "FE"
num_steps = [520*2**i for i in range(4)]
num_elements = [40*2**i for i in range(4)]

table(num_steps = num_steps, num_elements = num_elements, advance_method = advance_method,
    two_d = two_d, mono = mono)
