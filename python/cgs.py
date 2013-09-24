import numpy as np
import pandas as pd
import qgi
import sys
import yaml

N = qgi.N

default_params = dict(
    dh = 1e-3,
    graphs = [],
    s = [],
    h = [],
    cg = dict(
        eps = 1e-12,
        max_iter = 300
        )
    )

param_file = sys.argv[1] if len(sys.argv) > 1 else 'params.yaml'
stream = file(param_file, 'r')
p = default_params
p.update(yaml.load(stream))
stream.close()

delta = np.empty(qgi.D, dtype=np.dtype('d'))
resid = np.empty(qgi.D, dtype=np.dtype('d'))
psi = np.empty(qgi.D, dtype=np.dtype('d'))

def minimize(s, d, psi):
    return qgi.minimize_energy(s, d, p['cg']['eps'], p['cg']['max_iter'],
            psi, delta, resid)

data = []

for igraph, graph in enumerate(p['graphs']):

    graph = list('{:0>{width}b}'.format(int(graph, 16), width=N*(N-1)/2))
    graph = np.array(graph).astype(np.dtype('i'))
    psi0 = np.random.normal(size=qgi.D)
    h0prev = 0.

    for h0 in p['h']:
        for s in p['s']:

            h = h0 * np.ones(qgi.N)
            d = qgi.compute_diagonals(graph)
            qgi.update_diagonals(h0 - h0prev, d)
            energy, psi0 = minimize(s, d, psi0)
            mz = qgi.mag_z(psi0)
            mx = qgi.mag_x(psi0)
            q2 = qgi.overlap(psi0)
            
            q2p = 0.
            
            for j in range(qgi.N):
                qgi.update_diagonals_1(j, 0.5 * p['dh'], d)
                np.copyto(psi, psi0)
                _, psi = minimize(s, d, psi)
                fj = qgi.sigma_z(psi)
                
                qgi.update_diagonals_1(j, -p['dh'], d)
                np.copyto(psi, psi0)
                _, psi = minimize(s, d, psi)
                fj = (fj - qgi.sigma_z(psi)) / p['dh']
                q2p += np.sum(fj**2)
                
                qgi.update_diagonals_1(j, 0.5 * p['dh'], d)
                
            q2p = np.sqrt(q2p) / qgi.N
            data.append((igraph, h0, s, energy, mz, mx, q2, q2p))
            
data = pd.DataFrame.from_records(data,
    columns=['igraph', 'h', 's', 'energy', 'mz', 'mx', 'q2', 'q2p'])

params = [qgi.__version__, N, p['dh'], p['graphs']]
params = pd.Series(params, ['version', 'N', 'dh', 'graphs'])
data.set_index(['igraph', 'h', 'dh', 's'], inplace=True)

params.to_hdf(p['output_file'], 'params')
data.to_hdf(p['output_file'], 'data')
