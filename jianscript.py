import pandas
import numpy as np
import sys

filename = sys.argv[1] 
d = pandas.read_csv(filename, sep=',')
cols = [i for i in d.columns if i[0]=="N"]
new_cols = [i.replace("N","").replace("_","->").upper() for i in cols]

R1 = d.loc[d['read']==1, cols].sum(axis=0)

results = np.hstack([R1[i:i+4].values / R1[i:i+4].sum() for i in [0,4,8,12]])

for i,j in zip(new_cols, results):
    print(i,":",j)

