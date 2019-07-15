#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/NetBeansProjects/inapik/inapi0d0.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)

t=np.arange(-70,-10,0.1)
vnc=(0-8*(t+80)-20*(t-60)/(1+np.exp((-20-t)/15)))/(9*(t+90))
nnc=1/(1+np.exp((-25-t)/5))

nvalues=21
vvalues=21
points=51
file2=open('/home/richard/outhome/phaseliness0.txt',"r")
x,y=[],[]
for k in file2:
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

for vi in range(0,vvalues):
	for ni in range(0,nvalues):
		plt.plot(yv[vi][ni][:],xv[vi][ni][:])

plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
plt.plot(t,vnc,label='v-nullcline')
plt.plot(t,nnc,label='n-nullcline')
for vi in range(0,vvalues):
	for ni in range(0,nvalues):
		plt.plot(yv[vi][ni][:],xv[vi][ni][:])
#plt.gcf().subplots_adjust(left=0.15)
plt.plot(ax,ay)
plt.legend()
plt.savefig('inapburstwnpp.pdf')
