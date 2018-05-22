from conv_rates import *

two_d = 1
mono = True
advance_method = "OS"
num_steps = [10*4**i for i in range(3)]
num_elements = [60*2**i for i in range(3)]

table(num_steps = num_steps, num_elements = num_elements,
    advance_method = "OS", two_d = two_d, mono = mono)
