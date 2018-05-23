from conv_rates import *
from efficiency import error_

two_d = 0
mono = True
advance_method = "OS"
num_steps = [40*2**i for i in range(4)]
num_elements = [520*2**i for i in range(4)]

table(num_steps = num_steps, num_elements = num_elements,
    advance_method = advance_method, two_d = two_d, mono = mono)
