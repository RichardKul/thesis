#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

eqreltime=np.zeros((4,20))
breltime=np.zeros((4,20))
ii=0
for x in [50]:
	times,avalues,j2values,state=[],[],[],[]
	for y in range(15,16):
		file=open('/home/richard/mainp2test/time15atest%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			times.append(float(row[0]))
			avalues.append(float(row[1]))
			j2values.append(float(row[2]))
			state.append(float(row[3]))
		timesa=np.array(times)
		statea=np.array(state)
		countb=0
		counteq=0
		if statea[0]<0.5:
			intb=np.zeros(4999)
			inteq=np.zeros(5000)
			for z in range(0,5000):
				if timesa[2*z+1]>timesa[2*z]:
					inteq[z]=timesa[2*z+1]-timesa[2*z]
					counteq=counteq+1
				else:
					break	
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]				
			for z in range(0,4999):
				if timesa[2*z+2]>timesa[2*z+1]:
					intb[z]=timesa[2*z+2]-timesa[2*z+1]
					countb=countb+1
				else:
					break
			intbt=np.zeros(countb)
			for c in range(0,countb):
				intbt[c]=intb[c]
		else:
			inteq=np.zeros(4999)
			intb=np.zeros(5000)
			for z in range(0,5000):
				if timesa[2*z+1]>timesa[2*z]:
					intb[z]=timesa[2*z+1]-timesa[2*z]
					countb=countb+1
				else:
					break	
			intbt=np.zeros(countb)
			for c in range(0,countb):
				intbt[c]=intb[c]				
			for z in range(0,4999):
				if timesa[2*z+2]>timesa[2*z+1]:
					inteq[z]=timesa[2*z+2]-timesa[2*z+1]
					counteq=counteq+1
				else:
					break
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]	
				
		plt.hist(intbt, bins=50)
		plt.title("distribution of bursting time intervals")
		plt.savefig('bdist%d%d.pdf' %(x,y))
		plt.figure()
		plt.hist(inteqt, bins=50)
		plt.title("distribution of equilibrium time intervals")
		plt.savefig('eqdist%d%d.pdf' %(x,y))
		#eqtot=np.sum(inteq)
		#btot=np.sum(intb)
		#eqrel=eqtot/(eqtot+btot)
		#brel=btot/(eqtot+btot)
		#eqreltime[ii][y]=eqrel
		#breltime[ii][y]=brel
	ii=ii+1	

#xs=np.arange(-0.75,4.25,0.25)
#plt.xlabel('bias current I')
#plt.ylabel('proportion of simulation time')
#plt.plot(xs,eqreltime[0,:],label='D=1,burst')
#plt.plot(xs,eqreltime[1,:],label='D=2,burst')
#plt.plot(xs,eqreltime[2,:],label='D=3,burst')
#plt.plot(xs,eqreltime[3,:],label='D=4,burst')
#plt.plot(xs,breltime[0,:],label='D=1,equilibrium')
#plt.plot(xs,breltime[1,:],label='D=2,equilibrium')
#plt.plot(xs,breltime[2,:],label='D=3,equilibrium')
#plt.plot(xs,breltime[3,:],label='D=4,equilibrium')
#plt.legend()
#plt.savefig('inapiktimes.pdf')
