# Template parameters file

p = get_params()

p.dh = 1e-3

# G_13 and G_13^' from the paper by Vinci et al.
p.graphs = [
    '98401a4024000c14400',
    '000010102c1c0068c21'
]

p.s = np.linspace(0.02, 0.98, 49)
p.h = np.logspace(-2, 1, 52)
