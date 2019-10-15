#!/bin/python

import sys, re

def colored(character):
    return '\x1b[1;34;47m' + character + '\033[0m'

F1, F2 = open(sys.argv[1]), open(sys.argv[2])
f1,f2 = [i.strip('\n') for i in F1], [i.strip('\n') for i in F2]

n=0
firstline, secondline = [],[]
for i in range(1, len(f1[1])+1):
    if not i%10: 
        n+=1
        firstline.append(colored(str(n)))
        secondline.append(colored('0'))
    else:
        firstline.append(" ")
        secondline.append(str(i%10))

print(''.join(firstline))
print(''.join(secondline))
print(''.join(["-"]*len(f1[1])))

for w1,w2 in zip(f1,f2):
    pc = []

    for i,j in zip(w1,w2):
       if i==j:
           pc.append(i)
       else:
           pc.append(colored(i))

    print('{:76} | {}'.format(''.join(pc),w2))
