#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/NetBeansProjects/simplemodel/simplemodel5d2.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)

plt.figure()
plt.ylabel('gating variable')
plt.xlabel('membrane voltage V')
plt.plot(ax,ay)
plt.savefig('simplemodel5d2p.pdf')