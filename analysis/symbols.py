symbols = {'energy': 'E_0',
           'mz'    : 'm_z',
           'mx'    : 'm_x',
           'q2'    : 'Q_2'}

def symbol(s): return symbols[s] if s in symbols else s
