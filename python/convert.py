import qgi
import sys

N = qgi.N

for f in sys.argv[1:]:
    A = qgi.read_amatrix(f)
    b = ''.join((str(c) for c in qgi.tobitstr(A)))
    print '%0*x' % (N*(N-1)/8, int(b, 2))

