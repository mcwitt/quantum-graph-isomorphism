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
    df = pandas.read_table(files[0], sep=' *', index_col=['graph', 'h', 's'])
    for f in files[1:]:
        df = pandas.concat((df, pandas.read_table(f, sep=' *',
            index_col=['graph', 'h', 's'])))
    graphs = df.index.get_level_values('graph').unique()
    return graphs, df

def plot_heatmap(df, logy=True):
    x, y = np.meshgrid(df.columns, df.index)
    plt.pcolormesh(x, y, df.values)
    #plt.pcolormesh(x, y, df, norm=LogNorm())
    if logy: plt.yscale('log')
    plt.colorbar()
