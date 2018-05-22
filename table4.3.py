from conv_rates_monodomene import *

two_d = 0
mono = True
advance_method = "FE"
num_steps = [520*2**i for i in range(4)]
num_elements = [40*2**i for i in range(4)]

table(num_steps = num_steps, num_elements = num_elements, advance_method = "FE",
    two_d = False, mono = mono)
