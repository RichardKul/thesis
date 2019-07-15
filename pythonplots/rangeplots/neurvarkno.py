#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

vec=np.zeros(20)

col,colx=[],[]	
for y in range(1,21):
	file=open('/home/richard/outhome/dkno%d.txt' % (y),"r")
	for k in file:
		row=k.split()
		col.append(float(row[1]))
cola=np.array(col)
minv=np.amin(cola)
maxv=np.amax(cola)
x=np.arange(0.5,10.5,0.5)
#fig=plt.figure()
plt.yscale('log')
plt.plot(x,cola)
plt.yticks(np.arange(2450,2800,300))
plt.xlabel('time constant kn')
plt.ylabel('$D_{eff}$')
plt.savefig('dneurknovar.pdf')


