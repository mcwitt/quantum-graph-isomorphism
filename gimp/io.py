import numpy as np

def load_graph(filename):
    f = open(filename, 'r')
    a = []
    for line in f: a.append(list(line.strip()))
    a = np.array(a, dtype=int)
    assert(a.shape[0] == a.shape[1])
    return a

