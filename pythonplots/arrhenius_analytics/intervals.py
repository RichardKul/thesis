#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

eqreltime=np.zeros((4,20))
breltime=np.zeros((4,20))
ii=0
for x in [10,20,30,40]:
	times,state=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/times11a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			times.append(float(row[0]))
			state.append(float(row[1]))
		timesa=np.array(times)
		statea=np.array(state)
		if statea[0]<0.5:
			intb=np.zeros(4999)
			inteq=np.zeros(5000)
			for z in range(0,5000):
				inteq[z]=timesa[2*z+1]-timesa[2*z]
			for z in range(0,4999):
				intb[z]=timesa[2*z+2]-timesa[2*z+1]
		else:
			inteq=np.zeros(4999)
			intb=np.zeros(5000)
			for z in range(0,5000):
				intb[z]=timesa[2*z+1]-timesa[2*z]
			for z in range(0,4999):
				inteq[z]=timesa[2*z+2]-timesa[2*z+1]
		plt.hist(intb, bins=50)
		plt.title("distribution of bursting time intervals")
		plt.savefig('bdist%d%d.pdf' %(x,y))
		plt.hist(inteq, bins=50)
		plt.title("distribution of equilibrium time intervals")
		plt.savefig('eqdist%d%d.pdf' %(x,y))
		eqtot=np.sum(inteq)
		btot=np.sum(intb)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqreltime[ii][y]=eqrel
		breltime[ii][y]=brel
	ii=ii+1	

xs=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current I')
plt.ylabel('proportion of simulation time')
plt.plot(xs,eqreltime[0,:],label='D=1,burst')
plt.plot(xs,eqreltime[1,:],label='D=2,burst')
plt.plot(xs,eqreltime[2,:],label='D=3,burst')
plt.plot(xs,eqreltime[3,:],label='D=4,burst')
plt.plot(xs,breltime[0,:],label='D=1,equilibrium')
plt.plot(xs,breltime[1,:],label='D=2,equilibrium')
plt.plot(xs,breltime[2,:],label='D=3,equilibrium')
plt.plot(xs,breltime[3,:],label='D=4,equilibrium')
plt.legend()
plt.savefig('inapiktimes.pdf')
