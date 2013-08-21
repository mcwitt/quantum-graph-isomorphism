import pandas
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

symbols = {'energy': 'E_0',
           'mz'    : 'm_z',
           'mx'    : 'm_x',
           'q2'    : 'Q_2',
           'q2p'   : 'Q_2^{\prime}',
           'res2'  : 'r^2',
           'iterations': r'\mathrm{iterations}'}

def symbol(s): return symbols[s] if s in symbols else s

def load_files(files):
    index = ['ver', 'graph', 'h', 'dh', 's']
    df = pandas.read_table(files[0], sep=' *', index_col=index)
    for f in files[1:]:
        df = pandas.concat((df,
            pandas.read_table(f, sep=' *', index_col=index)))
    return df

def plot_heatmap(df, logy=True):
    x, y = np.meshgrid(df.columns, df.index)
    plt.pcolormesh(x, y, df.values)
    #plt.pcolormesh(x, y, df.values, norm=LogNorm())
    if logy: plt.yscale('log')
    plt.colorbar()
