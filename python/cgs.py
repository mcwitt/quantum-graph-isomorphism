import numpy as np
import pandas as pd
import qgi
import sys
from collections import namedtuple

N = qgi.N

def get_params():
    p = namedtuple('params', 'graphs s h dh cg output_file')
    p.cg = namedtuple('cg', 'max_iter eps')

    p.dh = 1e-3
    p.graphs = []
    p.s = []
    p.h = []
    p.cg.eps = 1e-12
    p.cg.max_iter = 300

    return p

param_file = sys.argv[1] if len(sys.argv) > 1 else 'params.py'
execfile(param_file)    # load parameters

delta = np.empty(qgi.D, dtype=np.dtype('d'))
resid = np.empty(qgi.D, dtype=np.dtype('d'))
psi = np.empty(qgi.D, dtype=np.dtype('d'))
data = []

for igraph, graph in enumerate(p.graphs):

    graph = list('{:0>{width}b}'.format(int(graph, 16), width=N*(N-1)/2))
    graph = np.array(graph).astype(np.dtype('i'))
    psi0 = np.random.normal(size=qgi.D)
    h0prev = 0.

    for h0 in p.h:
        for s in p.s:
            h = h0 * np.ones(qgi.N)
            d = qgi.compute_diagonals(graph)
            qgi.update_diagonals(h0 - h0prev, d)

            energy, psi0 = qgi.minimize_energy(s, d, p.cg.eps, p.cg.max_iter,
                    psi0, delta, resid)

            mz = qgi.mag_z(psi0)
            mx = qgi.mag_x(psi0)
            q2 = qgi.overlap(psi0)
            
            q2p = 0.
            
            for j in range(qgi.N):
                qgi.update_diagonals_1(j, 0.5 * p.dh, d)
                np.copyto(psi, psi0)
                _, psi = qgi.minimize_energy(s, d, p.cg.eps, p.cg.max_iter,
                        psi, delta, resid)
                fj = qgi.sigma_z(psi)
                
                qgi.update_diagonals_1(j, -p.dh, d)
                np.copyto(psi, psi0)
                _, psi = qgi.minimize_energy(s, d, p.cg.eps, p.cg.max_iter,
                        psi, delta, resid)
                fj = (fj - qgi.sigma_z(psi)) / p.dh
                q2p += np.sum(fj**2)
                
                qgi.update_diagonals_1(j, 0.5 * p.dh, d)
                
            q2p = np.sqrt(q2p) / qgi.N
            data.append((igraph, h0, s, energy, mz, mx, q2, q2p))
            
data = pd.DataFrame.from_records(data,
    columns=['igraph', 'h', 's', 'energy', 'mz', 'mx', 'q2', 'q2p'])

params = [qgi.__version__, p.dh, p.graphs]
params = pd.Series(params, ['version', 'dh', 'graphs'])
data.set_index(['igraph', 'h', 'dh', 's'], inplace=True)

params.to_hdf(p.output_file, 'params')
data.to_hdf(p.output_file, 'data')
