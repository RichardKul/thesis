#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


date='firingratestimesh'
D=[30,35,40,45,50]
dvalues=len(D)
istart=-0.05
irange=0.5
ivalues=50
istep=irange/ivalues
vec=np.zeros((dvalues,ivalues))
ii=0

file=open('/home/richard/outhome/%s.txt' %date,"r")
for k in file:
	col=[]
	row=k.split()
	for j in range(0,dvalues):
		col.append(float(row[j]))
	cola=np.array(col)	
	for z in range(0,dvalues):
		vec[z][ii]=cola[z]
	ii=ii+1

ivar=np.arange(istart,istart+irange,istep)
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.yscale('log')
#plt.xscale('log')
for n in range(0,dvalues):
	plt.plot(ivar,vec[n,:],label='D=%s' %D[n])

plt.savefig('firingratesh.pdf')


