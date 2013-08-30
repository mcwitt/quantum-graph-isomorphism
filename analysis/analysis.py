import pandas
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

symbols = {'energy': 'E_0',
           'mz'    : 'm_z',
           'mx'    : 'm_x',
           'q2'    : 'Q_2',
           'q2p'   : 'Q_2^{\prime}',
           'resid' : '|r|',
           'iterations': r'\mathrm{iterations}'}

measurements = ['q2p', 'q2', 'mz', 'mx', 'energy']

def symbol(s): return symbols[s] if s in symbols else s

def read_output(files):
    kwargs = dict(sep='[ |,] *', index_col=['ver', 'graph', 'h', 'dh', 's'])
    df = pandas.read_table(files[0], **kwargs)
    for f in files[1:]:
        df = pandas.concat((df, pandas.read_table(f, **kwargs)))
    return df

def plot_heatmap(df, logy=True, xylab=['$s$', '$h$']):
    x, y = np.meshgrid(df.columns, df.index)
    plt.pcolormesh(x, y, df.values)
    #plt.pcolormesh(x, y, df.values, norm=LogNorm())
    if logy: plt.yscale('log')
    plt.xlabel(xylab[0])
    plt.ylabel(xylab[1])
    plt.colorbar()
