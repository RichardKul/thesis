#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

vec=np.zeros((7,20))
ii=0

for x in [27,28,29,30,31,32,33]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	ii=ii+1

xs=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current I')
plt.ylabel('Fano factor')
#plt.yscale('log')
plt.xlim(0.5, 1.5)
plt.ylim(400, 700)
plt.plot(xs,vec[0,:],label='D=2.7')
plt.plot(xs,vec[1,:],label='D=2.8')
plt.plot(xs,vec[2,:],label='D=2.9')
plt.plot(xs,vec[3,:],label='D=3')
plt.plot(xs,vec[4,:],label='D=3.1')
plt.plot(xs,vec[5,:],label='D=3.2')
plt.plot(xs,vec[6,:],label='D=3.3')

plt.legend()
plt.savefig('fneurpmaxclose.pdf')


