import numpy as np
import pandas as pd
import qgi
import re
import scipy.sparse.linalg as arp
import sys

graph_files = sys.argv[1:]
output_file = 'output/{}-arpack.out'
s_vals = np.linspace(0.01, 0.99, 50)
h0vals = np.logspace(-2, 1, 49)
dh = 1e-3

def diagonalize(H, v0):
    energy, psi = arp.eigsh(H, v0=v0, k=1, which='SA')
    return energy[0], psi.T[0]

graphs = [qgi.load_amatrix(g) for g in graph_files]
for g in graphs: assert len(g) == qgi.N
output = []

for graph, graph_file in zip(graphs, graph_files):
    psi0 = np.random.normal(size=qgi.D)
    for h0 in h0vals:
        for s in s_vals:
            h = h0 * np.ones(qgi.N)
            H = qgi.hamiltonian(graph, h, s)
            energy, psi0 = diagonalize(H, psi0)
            mz = qgi.mag_z(psi0)
            mx = qgi.mag_x(psi0)
            q2 = qgi.overlap(psi0)
            
            q2p = 0.
            
            for j in range(qgi.N):
                h[j] += 0.5 * dh
                H = qgi.hamiltonian(graph, h, s)
                _, psi = diagonalize(H, psi0)
                fj = qgi.sigma_z(psi)
                
                h[j] -= dh
                H = qgi.hamiltonian(graph, h, s)
                _, psi = diagonalize(H, psi0)
                fj = (fj - qgi.sigma_z(psi)) / dh
                q2p += np.sum(fj**2)
                
                h[j] += 0.5 * dh
                
            q2p = np.sqrt(q2p) / qgi.N
            graph_file = re.sub(r'.*/(.*)$', r'\1', graph_file)
            output.append((graph_file, h0, s, energy, mz, mx, q2, q2p))
            
df = pd.DataFrame.from_records(output,
    columns=['graph', 'h', 's', 'energy', 'mz', 'mx', 'q2', 'q2p'])

df['ver'] = 'ARPACK'
df['dh'] = dh
df = df.set_index(['ver', 'graph', 'h', 'dh', 's'])
graph_file = re.sub(r'.*/(.*)$', r'\1', graph_files[0])
df.to_csv(output_file.format(graph_file))
