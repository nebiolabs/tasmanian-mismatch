#!/bin/python

import sys,os
import pandas 
import pickle
import numpy as np


for n,i in enumerate(sys.argv):
    if i=='-o':
        output_filename = sys.argv[n+1]

if not 'output_filename' in globals():
    exit("should specify -o output_file_name")

data = {}
for line in sys.stdin:
    category, length = line.strip().split(' ')
    length = int(length)

    if length > 10000: 
        continue
    
    if category not in data:
        data[category] = []

    data[category].append(length)


bins = np.hstack([np.arange(0,500,20),10000])


histogram = {}
for k,v in data.items():
    hist, b_edges = np.histogram(v, bins=bins)
    histogram[k] = hist


df = pandas.DataFrame(histogram)
df['bins'] = bins[1:]
df = df.set_index('bins')

df.to_csv(output_filename)

#with open('test.pkl', 'wb') as f:
#    pickle.dump(histogram, f)

