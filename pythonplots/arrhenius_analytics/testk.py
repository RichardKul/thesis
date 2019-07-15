#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

vec=np.zeros((3,2))
ii=0

for x in [0.1,0.2,0.3]:
    file=open('/home/richard/paratest/test%.1f.txt' % (x),"r")
    col=[]
    for k in file:
        col.append(k)
    for y in range(0,2):
        cola=np.array(col)
        vec[ii][y]=cola[y]
    ii=ii+1
print(vec[0][0],vec[1][1],vec[2][0])

