#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

vec=np.zeros((2,20))
ii=0

col=[]	
for y in range(1,21):
	file=open('/home/richard/outhome/d20%d.txt' % (y),"r")
	for k in file:
		row=k.split()
		col.append(row[1])
cola=np.array(col)
for z in range(0,20):
	vec[ii][z]=cola[z]
ii=ii+1
col=[]	
for y in range(1,21):
	file=open('/home/richard/outhome/d20a20%d.txt' % (y),"r")
	for k in file:
		row=k.split()
		col.append(row[1])
cola=np.array(col)
for z in range(0,20):
	vec[ii][z]=cola[z]


xs=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
plt.yscale('log')
plt.plot(xs,vec[0,:],label='D=2,old')
plt.plot(xs,vec[1,:],label='D=2,new')

plt.legend()
plt.savefig('d2neurpcomp.pdf')


