import fileinput
import numpy as np

# read in adjacency matrix

a = []

for line in fileinput.input():
    a.append(list(line.strip()))

a = np.array(a, dtype=int)

# construct Hamiltonian matrix


