import numpy as np

def spin(i, j): return (((1 << j) & i) >> j)*2 - 1

def neighbor(i, j): return (1 << j) ^ i

if __name__=='__main__':

    ## read in adjacency matrix ##
    import fileinput

    a = []

    for line in fileinput.input():
        a.append(list(line.strip()))

    a = np.array(a, dtype=int)
    n = a.shape[0]
    assert(a.shape[1] == n)

    ## construct Hamiltonian matrix ##

    d = 2**n    # dimension of Hilbert space

    H = np.empty((d, d))

    for i in xrange(d):
        for j in xrange(i):
            pass






