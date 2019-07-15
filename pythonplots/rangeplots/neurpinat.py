#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

vec=np.zeros((4,20))
ii=0

for x in [15,20,30,40]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/inatoptd%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	ii=ii+1

xs=np.arange(-0.9,1.1,0.1)
plt.xlabel('bias current I')
plt.ylabel('$D_eff$')
#plt.yscale('log')
plt.plot(xs,vec[0,:],label='D=1.5')
plt.plot(xs,vec[1,:],label='D=2')
plt.plot(xs,vec[2,:],label='D=3')
plt.plot(xs,vec[3,:],label='D=4')
#plt.plot(xs,vec[4,:],label='D=4')
#plt.plot(xs,vec[5,:],label='D=5')
#plt.plot(xs,vec[6,:],label='D=6')

plt.legend()
plt.savefig('dneurinat.pdf')


