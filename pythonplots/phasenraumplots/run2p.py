#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/NetBeansProjects/inatbursting/detinati0noise2.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)

plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
plt.plot(ax,ay)
plt.savefig('detinati0noise2p.pdf')
