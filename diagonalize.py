import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh

def spin(i, j): return (((1 << j) & i) >> j)*2 - 1

def neighbor(i, j): return (1 << j) ^ i

def hamiltonian(a, h, s):
    N = a.shape[0]
    d = 2**N
    H = lil_matrix((d, d))

    # driver Hamiltonian #

    for i in xrange(d):
        for j in xrange(N):
            H[i, neighbor(i, j)] = 1. - s

    # problem Hamiltonian #

    for i in xrange(d):

        for n in xrange(N):

            s_n = spin(i, n)
            H[i, i] += s * h[n] * s_n

            for m in xrange(n):
                if a[n, m] == 1:
                    H[i, i] += s * s_n * spin(i, m)

    return H


if __name__=='__main__':

    ## read in adjacency matrix ##
    import fileinput

    a = []

    for line in fileinput.input():
        a.append(list(line.strip()))

    a = np.array(a, dtype=int)
    assert(a.shape[0] == a.shape[1])
    h = np.ones(a.shape[0])
    H = hamiltonian(a, h, s=0.9)
    print eigsh(H, k=1)[0]
