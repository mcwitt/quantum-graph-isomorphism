import numpy as np
import pandas as pd
import qgi

N = qgi.N

def run_cgs(graphs, range_s, range_h, dh=0.001, eps=1e-12, itermax=300):

    # allocate memory for arrays
    delta = np.empty(qgi.D, dtype=np.dtype('d'))
    resid = np.empty(qgi.D, dtype=np.dtype('d'))
    psi = np.empty(qgi.D, dtype=np.dtype('d'))

    def minimize(s, d, psi):
        return qgi.minimize_energy(s, d, eps, itermax, psi, delta, resid)

    output = []

    for igraph, graph in enumerate(graphs):

        graph = list('{:0>{width}b}'.format(int(graph, 16), width=N*(N-1)/2))
        graph = np.array(graph).astype(np.dtype('i'))
        psi0 = np.random.normal(size=qgi.D)
        h0prev = 0.

        for h0 in range_h:
            for s in range_s:
                h = h0 * np.ones(N)
                d = qgi.compute_diagonals(graph)
                qgi.update_diagonals(h0 - h0prev, d)
                energy, psi0 = minimize(s, d, psi0)
                mz = qgi.mag_z(psi0)
                mx = qgi.mag_x(psi0)
                q2 = qgi.overlap(psi0)
                q2p = 0.
                
                for j in range(N):
                    qgi.update_diagonals_1(j, 0.5 * dh, d)
                    np.copyto(psi, psi0)
                    _, psi = minimize(s, d, psi)
                    fj = qgi.sigma_z(psi)
                    
                    qgi.update_diagonals_1(j, -dh, d)
                    np.copyto(psi, psi0)
                    _, psi = minimize(s, d, psi)
                    fj = (fj - qgi.sigma_z(psi)) / dh 
                    q2p += np.sum(fj**2)
                    
                    qgi.update_diagonals_1(j, 0.5 * dh, d)
                    
                q2p = np.sqrt(q2p) / N
                output.append((igraph, h0, s, energy, mz, mx, q2, q2p))
                
    return pd.DataFrame.from_records(output,
        columns=['igraph', 'h', 's', 'energy', 'mz', 'mx', 'q2', 'q2p'])

if __name__ == '__main__':
    import sys
    import yaml

    default_params = dict(
            N = N,
            dh = 0.001,
            graphs = [],
            range = dict(
                s = [0.02, 0.98, 0.02],
                log_h = [-2, 1, 0.03]
                ),
            cg = dict(
                eps = 1e-12,
                itermax = 300
                ),
            output_file = 'out.h5'
            )

    param_file = sys.argv[1] if len(sys.argv) > 1 else 'params.yaml'
    stream = file(param_file, 'r')
    p = default_params
    p.update(yaml.load(stream))
    stream.close()

    if p['N'] != N:
        msg = '{} error: expected {} vertices, got {}'.format(sys.argv[0], N, p.N)
        sys.stderr.write(msg)
        sys.exit(1)

    range_s = np.arange(*p['range']['s'])
    range_h = 10**(np.arange(*p['range']['log_h']))
    output = run_cgs(p['graphs'], range_s, range_h, p['dh'], **p['cg'])
    output.set_index(['igraph', 'h', 'dh', 's'], inplace=True)
    params = [qgi.__version__, N, p['dh'], p['graphs']]
    params = pd.Series(params, ['version', 'N', 'dh', 'graphs'])

    # write data to HDF file
    params.to_hdf(p['output_file'], 'params')
    output.to_hdf(p['output_file'], 'output')
