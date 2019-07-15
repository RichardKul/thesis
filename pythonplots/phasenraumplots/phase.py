#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

nvalues=21
vvalues=21
points=51
file=open('/home/richard/phasenraum/phasebw1.txt',"r")
x,y=[],[]
for k in file:
    row=k.split()
    x.append(row[2])
    y.append(row[3])
ax=np.array(x)
ay=np.array(y)
xv=np.zeros((vvalues,nvalues,points))
yv=np.zeros((vvalues,nvalues,points))
for vi in range(0,vvalues):
	for ni in range(0,nvalues):
		for pi in range(0,points):
			xv[vi][ni][pi]=ax[(ni*vvalues+vi)*points+pi]	
			yv[vi][ni][pi]=ay[(ni*vvalues+vi)*points+pi]

plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
for vi in range(0,vvalues):
	for ni in range(0,nvalues):
		plt.plot(yv[vi][ni][:],xv[vi][ni][:])
		#plt.arrow(yv[vi][ni][points-2],xv[vi][ni][points-2],yv[vi][ni][points-1]-yv[vi][ni][points-2],xv[vi][ni][points-1]-xv[vi][ni][points-2],head_width=10*0.003,head_length=15*0.003)
plt.savefig('phasebw1.pdf')
