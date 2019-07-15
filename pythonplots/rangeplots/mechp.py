#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


col=[]	
for y in range(1,31):
	file=open('/home/richard/outhome/mechd%d.txt' % (y),"r")
	for k in file:
		row=k.split()
		col.append(row[1])
cola=np.array(col)
	
xs=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current I')
plt.ylabel('$D_eff$')
plt.yscale('log')
plt.plot(xs,vec[4],label='D=1')
plt.plot(xs,vec[0,:],label='D=1.2')
plt.plot(xs,vec[1,:],label='D=2')
#plt.plot(xs,vec[[2],:],label='D=2')
plt.plot(xs,vec[2,:],label='D=3')
plt.plot(xs,vec[3,:],label='D=4')
#plt.plot(xs,vec[[5],:],label='D=5')
#plt.plot(xs,vec[[6],:],label='D=6')

plt.legend()
plt.savefig('dneurp.pdf')


