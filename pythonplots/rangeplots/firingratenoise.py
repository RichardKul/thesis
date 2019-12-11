#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


date='firingratestimesh'
D=[150,200,300,500]
dvalues=len(D)
T=25000
istart=-21
irange=15
ivalues=51
istep=irange/ivalues
vec=np.zeros((dvalues,ivalues))
ii=0

file=open('/home/richard/NetBeansProjects/detrinzel/countrinzelnoise4.txt',"r")
for k in file:
	col=[]
	row=k.split()
	for j in range(0,dvalues+1):
		col.append(float(row[j]))
	cola=np.array(col)	
	for z in range(0,dvalues):
		vec[z][ii]=cola[z+1]
	ii=ii+1

ivar=np.arange(istart,istart+irange,istep)
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.yscale('log')
#plt.xscale('log')
for n in range(0,dvalues):
	plt.plot(ivar,vec[n,:]/T,label='D=%s' %(D[n]/10))
plt.legend()
plt.savefig('firingraterinzelnoise3.pdf')


