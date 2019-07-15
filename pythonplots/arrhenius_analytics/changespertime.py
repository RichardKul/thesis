#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def func(x, a, b):
	return a * np.exp(-b * x)

changespertime=np.zeros((5,20))
eqtottime=np.zeros((5,20))
eqreltime=np.zeros((5,20))
btottime=np.zeros((5,20))
breltime=np.zeros((5,20))
ii=0
dt=0.00001
N0=320000000
Neq0=50000000
Ndiff0=N0-Neq0
repetitions=500
runs0=10
for x in [12]:
	for y in [1,2,3,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20]:
		times,avalues,jrvalues,j2values,state=[],[],[],[],[]
		file=open('/home/richard/outhome/time26a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			avalues.append(float(row[1]))
			jrvalues.append(float(row[2]))
			j2values.append(float(row[3]))
			state.append(float(row[4]))
		j2valuesa=np.array(j2values)
		jrvaluesa=np.array(jrvalues)
		avaluesa=np.array(avalues)
		statea=np.array(state)
		countb=0
		counteq=0
		if statea[0]<0.5:
			intb=np.zeros(4999)
			inteq=np.zeros(5000)
			for z in range(0,5000):
				if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff0>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff0:
					inteq[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-jrvaluesa[2*z])*Ndiff0+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					counteq=counteq+1
				else:
					break	
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff0>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff0:
					intb[z]=(((avaluesa[2*z+2]-avaluesa[2*z+1])*repetitions+jrvaluesa[2*z+2]-jrvaluesa[2*z+1])*Ndiff0+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
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
				if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff0>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff0:
					intb[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-jrvaluesa[2*z])*Ndiff0+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					countb=countb+1
				else:
					break	
			intbt=np.zeros(countb)
			for c in range(0,countb):
				intbt[c]=intb[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff0>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff0:
					inteq[z]=(((avaluesa[2*z+2]-avaluesa[2*z+1])*repetitions+jrvaluesa[2*z+2]-jrvaluesa[2*z+1])*Ndiff0+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
					counteq=counteq+1
				else:
					break
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]	
		counts=counteq+countb
		countsrel=counts/Ndiff0
		changespertime[ii][y-1]=countsrel/(repetitions*runs0)
		#plt.figure()		
		#plt.hist(intbt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of bursting time intervals")
		#plt.savefig('bdistajrj2%d%d.pdf' %(x,y))
		#plt.figure()
		#plt.hist(inteqt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of equilibrium time intervals")
		#plt.savefig('eqdistajrj2%d%d.pdf' %(x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-1]=eqtot/counteq
		btottime[ii][y-1]=btot/countb
		eqreltime[ii][y-1]=eqrel
		breltime[ii][y-1]=brel
	changespertime[ii][3]=changespertime[ii][2]
	changespertime[ii][5]=changespertime[ii][4]
	btottime[ii][3]=btottime[ii][2]
	btottime[ii][5]=btottime[ii][4]
	eqtottime[ii][3]=eqtottime[ii][2]
	eqtottime[ii][5]=eqtottime[ii][4]
	breltime[ii][3]=breltime[ii][2]
	breltime[ii][5]=breltime[ii][4]
	eqreltime[ii][3]=eqreltime[ii][2]
	eqreltime[ii][5]=eqreltime[ii][4]
	ii=ii+1	
N=120000000
Neq=20000000
Ndiff=N-Neq
#repetitions=500
runs=10
for x in [20]:
	for y in range(1,21):
		times,avalues,jrvalues,j2values,state=[],[],[],[],[]
		file=open('/home/richard/outhome/time20a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			avalues.append(float(row[1]))
			jrvalues.append(float(row[2]))
			j2values.append(float(row[3]))
			state.append(float(row[4]))
		j2valuesa=np.array(j2values)
		jrvaluesa=np.array(jrvalues)
		avaluesa=np.array(avalues)
		statea=np.array(state)
		countb=0
		counteq=0
		if statea[0]<0.5:
			intb=np.zeros(4999)
			inteq=np.zeros(5000)
			for z in range(0,5000):
				if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff:
					inteq[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-jrvaluesa[2*z])*Ndiff+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					counteq=counteq+1
				else:
					break	
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff:
					intb[z]=(((avaluesa[2*z+2]-avaluesa[2*z+1])*repetitions+jrvaluesa[2*z+2]-jrvaluesa[2*z+1])*Ndiff+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
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
				if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff:
					intb[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-jrvaluesa[2*z])*Ndiff+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					countb=countb+1
				else:
					break	
			intbt=np.zeros(countb)
			for c in range(0,countb):
				intbt[c]=intb[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff:
					inteq[z]=(((avaluesa[2*z+2]-avaluesa[2*z+1])*repetitions+jrvaluesa[2*z+2]-jrvaluesa[2*z+1])*Ndiff+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
					counteq=counteq+1
				else:
					break
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]	
		
		counts=counteq+countb
		countsrel=counts/Ndiff
		changespertime[ii][y-1]=countsrel/(repetitions*runs)
		#plt.figure()		
		#plt.hist(intbt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of bursting time intervals")
		#plt.savefig('bdistajrj2%d%d.pdf' %(x,y))
		#plt.figure()
		#plt.hist(inteqt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of equilibrium time intervals")
		#plt.savefig('eqdistajrj2%d%d.pdf' %(x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-1]=eqtot/counteq
		btottime[ii][y-1]=btot/countb
		eqreltime[ii][y-1]=eqrel
		breltime[ii][y-1]=brel
	ii=ii+1	

N2=100000000
Neq2=10000000
Ndiff2=N2-Neq2
runs=50
repetitions=100
for x in [30,40,50]:
	for y in range(1,21):
		times,avalues,jrvalues,j2values,state=[],[],[],[],[]
		file=open('/home/richard/outhome/time1m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			avalues.append(float(row[1]))
			jrvalues.append(float(row[2]))
			j2values.append(float(row[3]))
			state.append(float(row[4]))
		j2valuesa=np.array(j2values)
		jrvaluesa=np.array(jrvalues)
		avaluesa=np.array(avalues)
		statea=np.array(state)
		countb=0
		counteq=0
		if statea[0]<0.5:
			intb=np.zeros(4999)
			inteq=np.zeros(5000)
			for z in range(0,5000):
				if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff2>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff2:
					inteq[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-jrvaluesa[2*z])*Ndiff2+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					counteq=counteq+1
				else:
					break	
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff2>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff2:
					intb[z]=(((avaluesa[2*z+2]-avaluesa[2*z+1])*repetitions+jrvaluesa[2*z+2]-jrvaluesa[2*z+1])*Ndiff2+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
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
				if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff2>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff2:
					intb[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-jrvaluesa[2*z])*Ndiff2+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
					countb=countb+1
				else:
					break	
			intbt=np.zeros(countb)
			for c in range(0,countb):
				intbt[c]=intb[c]				
			for z in range(0,4999):
				if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff2>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff2:
					inteq[z]=(((avaluesa[2*z+2]-avaluesa[2*z+1])*repetitions+jrvaluesa[2*z+2]-jrvaluesa[2*z+1])*Ndiff2+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
					counteq=counteq+1
				else:
					break
			inteqt=np.zeros(counteq)
			for c in range(0,counteq):
				inteqt[c]=inteq[c]	
		
		counts=counteq+countb
		countsrel=counts/Ndiff2
		if counts<9999:
			changespertime[ii][y-1]=countsrel/(repetitions*runs)
		else:
			changespertime[ii][y-1]=countsrel/(avaluesa[9999]*repetitions+jrvaluesa[9999]+j2valuesa[9999]/Ndiff2)					
		#plt.figure()		
		#plt.hist(intbt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of bursting time intervals")
		#plt.savefig('bdistajrj2%d%d.pdf' %(x,y))
		#plt.figure()
		#plt.hist(inteqt, bins=50)
		#plt.yscale('log')
		#plt.title("distribution of equilibrium time intervals")
		#plt.savefig('eqdistajrj2%d%d.pdf' %(x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-1]=eqtot/counteq
		btottime[ii][y-1]=btot/countb
		eqreltime[ii][y-1]=eqrel
		breltime[ii][y-1]=brel
	ii=ii+1
#N3=50000000
#Neq3=10000000
#Ndiff3=N3-Neq3
#for x in [30,40]:
#	for y in range(10,11):
#		times,avalues,j2values,state=[],[],[],[]
#		file=open('/home/richard/outhome/time17a%d%d.txt' % (x,y),"r")
#		for k in file:
#			row=k.split()
#			avalues.append(float(row[1]))
#			j2values.append(float(row[2]))
#			state.append(float(row[3]))
#		j2valuesa=np.array(j2values)
#		avaluesa=np.array(avalues)
#		statea=np.array(state)
#		countb=0
#		counteq=0
#		if statea[0]<0.5:
#			intb=np.zeros(4999)
#			inteq=np.zeros(5000)
#			for z in range(0,5000):
#				if avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff3>avaluesa[2*z]+j2valuesa[2*z]/Ndiff3:
#					inteq[z]=((avaluesa[2*z+1]-avaluesa[2*z])*Ndiff3+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
#					counteq=counteq+1
#				else:
#					break	
#			inteqt=np.zeros(counteq)
#			for c in range(0,counteq):
#				inteqt[c]=inteq[c]				
#			for z in range(0,4999):
#				if avaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff3>avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff3:
#					intb[z]=((avaluesa[2*z+2]-avaluesa[2*z+1])*Ndiff3+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
#					countb=countb+1
#				else:
#					break
#			intbt=np.zeros(countb)
#			for c in range(0,countb):
#				intbt[c]=intb[c]
#		else:
#			inteq=np.zeros(4999)
#			intb=np.zeros(5000)
#			for z in range(0,5000):
#				if avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff3>avaluesa[2*z]+j2valuesa[2*z]/Ndiff3:
#					intb[z]=((avaluesa[2*z+1]-avaluesa[2*z])*Ndiff3+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
#					countb=countb+1
#				else:
#					break	
#			intbt=np.zeros(countb)
#			for c in range(0,countb):
#				intbt[c]=intb[c]				
#			for z in range(0,4999):
#				if avaluesa[2*z+2]+j2valuesa[2*z+2]/Ndiff3>avaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff3:
#					inteq[z]=((avaluesa[2*z+2]-avaluesa[2*z+1])*Ndiff3+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
#					counteq=counteq+1
#				else:
#					break
#			inteqt=np.zeros(counteq)
#			for c in range(0,counteq):
#				inteqt[c]=inteq[c]	
#		
#		#plt.figure()		
#		#plt.hist(intbt, bins=50)
#		#plt.yscale('log')
#		#plt.title("distribution of bursting time intervals")
#		#plt.savefig('bdistaj2%d%d.pdf' %(x,y))
#		#plt.figure()
#		#plt.hist(inteqt, bins=50)
#		#plt.yscale('log')
#		#plt.title("distribution of equilibrium time intervals")
#		#plt.savefig('eqdistaj2%d%d.pdf' %(x,y))
#		eqtot=np.sum(inteqt)
#		btot=np.sum(intbt)
#		eqrel=eqtot/(eqtot+btot)
#		brel=btot/(eqtot+btot)
#		eqreltime[ii][0]=eqrel
#		breltime[ii][0]=brel
#	ii=ii+1	
#
#xs=np.arange(0.25,4.25,0.25)
#for k in range(0,20):
#	plt.figure()
#	xs=[1/2,1/3,1/4,1/5]
#	plt.xlabel('inverse noise intensity 1/D')
#	plt.ylabel('transition rate')
#	plt.yscale('log')
#	plt.plot(xs,1/breltime[:,k],label='burst to eq')
#plt.plot(xs,breltime[1,:],label='D=3,burst')
#plt.plot(xs,breltime[2,:],label='D=4,burst')
#plt.plot(xs,eqreltime[3,:],label='D=4,burst')
#	plt.plot(xs,1/eqreltime[:,k],label='eq to burst')
#plt.plot(xs,eqreltime[1,:],label='D=3,equilibrium')
#plt.plot(xs,eqreltime[2,:],label='D=4,equilibrium')
#plt.plot(xs,breltime[3,:],label='D=4,equilibrium')
#	plt.legend()
#	plt.savefig('arrhenius%d.pdf' %(k))
#for k2 in range(5,16):
#	plt.figure()
#	xs=[1/2,1/3,1/4,1/5]
#	plt.xlabel('inverse noise intensity 1/D')
#	plt.ylabel('transition rate')
#	plt.yscale('log')
#	plt.plot(xs,1/btottime[:,k2],'bo',label='burst to eq')
#plt.plot(xs,breltime[1,:],label='D=3,burst')
#plt.plot(xs,breltime[2,:],label='D=4,burst')
#plt.plot(xs,eqreltime[3,:],label='D=4,burst')
#	plt.plot(xs,1/eqtottime[:,k2],'ro',label='eq to burst')
#plt.plot(xs,eqreltime[1,:],label='D=3,equilibrium')
#plt.plot(xs,eqreltime[2,:],label='D=4,equilibrium')
#plt.plot(xs,breltime[3,:],label='D=4,equilibrium')
#	popt,pcov = curve_fit(func, xs, 1/btottime[:,k2])
#	plt.plot(np.array(xs), func(np.array(xs), *popt), 'b-',label='fit burst to eq: r_0=%5.3f, U_+=%5.3f' % tuple(popt))
#	params[0][k2-5]=popt[0]
#	params[1][k2-5]=popt[1]
#	popt,pcov = curve_fit(func, xs, 1/eqtottime[:,k2])
#	plt.plot(np.array(xs), func(np.array(xs), *popt), 'r-',label='fiteq to burst: r_0=%5.3f, U_-=%5.3f' % tuple(popt))
#	params[2][k2-5]=popt[0]
#	params[3][k2-5]=popt[1]
#	plt.legend()
#	plt.savefig('arrheniustot1mfit%d.pdf' %(k2))
#	eqfile = open('param%d.txt' % (k2),'w')
#	for k3 in range(0,4): 
#		eqfile.write('%.6f\n'%params[k3][k2-5]) 
#	eqfile.close() 
#	plt.figure()
#	plt.xlabel('inverse noise intensity 1/D')
#	plt.ylabel('correlation time')
#	plt.yscale('log')
#	plt.plot(xs,1/(1/btottime[:,k2]+1/eqtottime[:,k2]))
#	plt.savefig('cortime%d.pdf' %(k2))
#plt.figure()
#xold=np.arange(-0.75,4.25,0.25)
#plt.xlabel('bias current I')
#plt.ylabel('correlation time')
#plt.yscale('log')
#plt.plot(xold,1/(1/btottime[0,:]+1/eqtottime[0,:]),label='D=2')
#plt.plot(xold,1/(1/btottime[1,:]+1/eqtottime[1,:]),label='D=3')
#plt.plot(xold,1/(1/btottime[2,:]+1/eqtottime[2,:]),label='D=4')
#plt.plot(xold,1/(1/btottime[3,:]+1/eqtottime[3,:]),label='D=5')
#plt.savefig('altcortime.pdf')

#plt.figure()
#xnew=np.arange(0.25,3,0.25)
#plt.xlabel('bias current')
#plt.ylabel('prefactor')
#plt.plot(xnew,params[0,:],label='burst to eq')
#plt.plot(xnew,params[2,:],label='eq to burst')
#plt.legend()
#plt.savefig('prefac.pdf')
#plt.figure()
#plt.xlabel('bias current')
#plt.ylabel('potential barrier')
#plt.plot(xnew,params[1,:],label='burst to eq')
#plt.plot(xnew,params[3,:],label='eq to burst')
#plt.legend()
#plt.savefig('barrier.pdf')
#plt.figure()
xs1=np.arange(-0.75,4.25,0.25)
xs2=np.arange(0.2,4.2,0.2)	
plt.xlabel('bias current')
plt.ylabel('changes over time steps')
plt.yscale('log')
plt.plot(xs2,changespertime[0,:],label='D=1.2')
plt.plot(xs1,changespertime[1,:],label='D=2')
plt.plot(xs1,changespertime[2,:],label='D=3')
plt.plot(xs1,changespertime[3,:],label='D=4')
plt.plot(xs1,changespertime[4,:],label='D=5')
plt.legend()
plt.savefig('changespertime1m.pdf')
