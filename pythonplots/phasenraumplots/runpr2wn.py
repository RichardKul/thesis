#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

nvalues=21
vvalues=21
points=51
file=open('/home/richard/outhome/phaseliness6.txt',"r")
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

t=np.arange(-80,0,0.01)
vnc=(6-8*(t+80)-20*(t-60)/(1+np.exp((-20-t)/15)))/(9*(t+90))
nnc=1/(1+np.exp((-25-t)/5))
plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
for vi in range(0,vvalues):
	for ni in range(0,nvalues):
		plt.plot(yv[vi][ni][:],xv[vi][ni][:])
		#plt.arrow(yv[vi][ni][points-2],xv[vi][ni][points-2],yv[vi][ni][points-1]-yv[vi][ni][points-2],xv[vi][ni][points-1]-xv[vi][ni][points-2],head_width=10*0.003,head_length=15*0.003)
plt.plot(t,vnc,color='black')
plt.plot(t,nnc,color='black')
plt.savefig('phaseliness6wn.pdf')

