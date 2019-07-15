#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/NetBeansProjects/inatbursting/detinatopt2i0d0.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax=np.array(x)
ay=np.array(y)

plt.figure()
plt.xlabel('time')
plt.ylabel('membrane voltage')
plt.plot(ax,ay)
plt.savefig('detinatopt2i0d0.pdf')
