from conv_rates import *

two_d = 1
mono = True
advance_method = "FE"
num_steps = [200*2**i for i in range(3)]
num_elements = [35*2**i for i in range(3)]

table(num_steps = num_steps, num_elements = num_elements,
    advance_method = advance_method, two_d = two_d, mono = mono)
