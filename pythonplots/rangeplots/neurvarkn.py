#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

vec=np.zeros(20)

col,colx=[],[]	
for y in range(1,21):
	file=open('/home/richard/outhome/dknn%d.txt' % (y),"r")
	for k in file:
		row=k.split()
		col.append(float(row[1]))
cola=np.array(col)
x=np.arange(0.5,10.5,0.5)
plt.xlabel('time constant kn')
plt.ylabel('$D_{eff}$')
plt.yscale('log')
plt.plot(x,cola)

plt.savefig('dneurknnvar.pdf')


