import pandas
import numpy as np
import sys

filename = sys.argv[1] 
d = pandas.read_csv(filename, sep=',')
cols = [i for i in d.columns if i[0]=="N"]
new_cols = [i.replace("N","").replace("_","->").upper() for i in cols]

R1 = d.loc[d['read']==1, cols].sum(axis=0)

alts = np.hstack([R1[i:i+4].values for i in [0,4,8,12]])
refs = np.hstack([4*[R1[i:i+4].sum()] for i in [0,4,8,12]])

results = alts/refs

print("{}\t{}\t{}\t{}".format('mismatch', 'ref_count', 'alt_count','rate'))

for i,j,k,l in zip(new_cols, refs, alts, results):
    print("{}\t{}\t{}\t{:.6f}".format(i,j,k,l))

