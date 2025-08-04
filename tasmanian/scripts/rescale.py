import pandas as pd
import numpy as np
from sys import argv

d = pd.read_csv('rescaling_matrix.csv')

A = list(d.columns[2:6])
T = list(d.columns[6:10])
C = list(d.columns[10:14])
G = list(d.columns[14:18])
all=['Na_t', 'Na_c', 'Na_g', 'Nt_a',  'Nt_c', 'Nt_g', 'Nc_a', 'Nc_t', 'Nc_g', 'Ng_a', 'Ng_t', 'Ng_c']

#d.loc[:, 'A'] = d.loc[:, A].sum(axis=1)
#d.loc[:, 'T'] = d.loc[:, T].sum(axis=1)
#d.loc[:, 'C'] = d.loc[:, C].sum(axis=1)
#d.loc[:, 'G'] = d.loc[:, G].sum(axis=1)

df = d.copy().astype(float)


for i in [A,C,G,T]:
    for j in i:
        df.loc[:,j] = df.loc[:,j].values / (df.loc[:, i].sum(axis=1).values + 1 ) # add 1 to avoid dividing by zero during tests.


df.loc[:, all] = 1 / 10**(df.loc[:, all])

df.loc[:, ['read', 'position'] + all].to_csv('rescaling_factors.csv', index=False)
