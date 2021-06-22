#!/usr/bin/env python

#A	C	A>C	184	0	0
#A	G	A>G	184	0	0
#A	T	A>T	184	0	0
#C	A	C>A	118	0	0
#C	G	C>G	118	0	0
#C	T	C>T	118	2	0.016667

import sys
import pandas
import matplotlib.pyplot as plt

d = []
for line in sys.stdin:
    ref, alt, subst, ref_count, alt_count, rate = line.strip().split("\t")
    d.append(line.strip().split("\t"))

d = pandas.DataFrame(d)
d.columns = ['ref', 'alt', 'subst', 'ref_count', 'alt_count', 'rate']
d.set_index('subst', inplace=True)

plt.bar(list(d.index), d.loc[:,'rate'])
plt.savefig('test.png', dpi=300)
