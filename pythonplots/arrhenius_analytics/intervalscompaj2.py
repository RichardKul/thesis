#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

eqreltime=np.zeros((1,20))
breltime=np.zeros((1,20))
ii=0
dt=0.00001
N=50000000
Neq=10000000
Ndiff=N-Neq
for x in [30]:
	for y in range(1,21):
		times,avalues,j2values,state=[],[],[],[]
		file=open('/home/richard/outhome/time17a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			avalues.append(float(row[1]))
			j2values.append(float(row[2]))
			state.append(float(row[3]))
		j2valuesa=np.array(j2values)
		avaluesa=np.array(avalues)
		statea=np.array(state)
		countb=0
		counteq=0
		if statea[0]<0.5:
			intb=np.zeros(4999)
			inteq=np.zeros(5000)
			for z in range(0,5000):
				if avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff>avaluesa[2*z]+j2valuesa[2*z]/Ndiff:
					inteq[z]=((avaluesa[2*z+1]-avaluesa[2*z])*Ndiff+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					counteq=counteq+1
				else:
					break	
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff>avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff:
					intb[z]=((avaluesa[2*z+2]-avaluesa[2*z+1])*Ndiff+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
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
				if avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff>avaluesa[2*z]+j2valuesa[2*z]/Ndiff:
					intb[z]=((avaluesa[2*z+1]-avaluesa[2*z])*Ndiff+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					countb=countb+1
				else:
					break	
			intbt=np.zeros(countb)
			for c in range(0,countb):
				intbt[c]=intb[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff>avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff:
					inteq[z]=((avaluesa[2*z+2]-avaluesa[2*z+1])*Ndiff+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
					counteq=counteq+1
				else:
					break
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]	
		
		#plt.figure()		
		#plt.hist(intbt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of bursting time intervals")
		#plt.savefig('bdistaj2%d%d.pdf' %(x,y))
		#plt.figure()
		#plt.hist(inteqt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of equilibrium time intervals")
		#plt.savefig('eqdistaj2%d%d.pdf' %(x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqreltime[ii][y-1]=eqrel
		breltime[ii][y-1]=brel
	ii=ii+1	

plt.figure()
xs=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current I')
plt.ylabel('proportion of simulation time')
plt.plot(xs,eqreltime[0,:],label='D=4,burst')
#plt.plot(xs,eqreltime[1,:],label='D=2,burst')
#plt.plot(xs,eqreltime[2,:],label='D=3,burst')
#plt.plot(xs,eqreltime[3,:],label='D=4,burst')
plt.plot(xs,breltime[0,:],label='D=4,equilibrium')
#plt.plot(xs,breltime[1,:],label='D=2,equilibrium')
#plt.plot(xs,breltime[2,:],label='D=3,equilibrium')
#plt.plot(xs,breltime[3,:],label='D=4,equilibrium')
plt.legend()
plt.savefig('inapiktimes173.pdf')
