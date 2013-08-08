symbols = {'energy': 'E_0',
           'mz'    : 'm_z',
           'mx'    : 'm_x',
           'q2'    : 'Q_2',
           'res2'  : 'r^2',
           'iterations': r'\mathrm{iterations}'}

def symbol(s): return symbols[s] if s in symbols else s
