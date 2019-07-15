#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

#vec=np.zeros((1,20))
#ii=0
#vec=np.zeros((1,7))

#for x in [8,10,12,15,20,30,40]:
for x in [3]:
	col,colx=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/drealtest30j22%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
#	for z in range(0,20):
#		vec[ii][z]=cola[z]
#	ii=ii+1

#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('$D_{eff}$')
plt.yscale('log')
#plt.xscale('log')
#plt.plot(xs,vec[0,:],label='D=0.8')
#plt.plot(xs,vec[1,:],label='D=1')
#plt.plot(xs,vec[2,:],label='D=1.2')
#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
plt.plot(colxa,cola,label='D=3e-3')

plt.legend()
plt.savefig('dneur30j223.pdf')


