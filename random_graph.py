import sys
import numpy as np
from numpy.random import random_integers

N = int(sys.argv[1])

a = np.empty((N, N))

for i in xrange(N):

    a[i, i] = 0

    for j in xrange(i):
        a[i, j] = a[j, i] = random_integers(0, 1)

for i in xrange(N):

    for j in xrange(N):
        sys.stdout.write('%d' % a[i, j])

    sys.stdout.write('\n')
