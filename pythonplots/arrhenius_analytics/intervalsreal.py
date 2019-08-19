#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def func(x, a, b):
	return a * np.exp(-b * x)
def func2(x, a, b):
	return a * x + b
def func3(x, a):
	return a
def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def deffana(t,r0p,r0m,ap,am,bp,bm):
	return (r(r0p,ap,bp,t,30)*r(r0m,am,bm,t,30))/((r(r0p,ap,bp,t,30)+r(r0m,am,bm,t,30))**3)

D=[20]
l=len(D)
Da=np.array(D)
d=Da[0]/100

#def deffana(t,r0p,r0m,ap,am,bp,bm):
#	return (r(r0p,ap,bp,t,0.5)*r(r0m,am,bm,t,0.5))/((r(r0p,ap,bp,t,0.5)+r(r0m,am,bm,t,0.5))**3) 

date='realfast11jjem2'
params2=np.zeros(6)
params=np.zeros((4,11))
changespertime=np.zeros((l,20))
eqtottime=np.zeros((l,20))
eqreltime=np.zeros((l,20))
btottime=np.zeros((l,20))
breltime=np.zeros((l,20))
countstot=np.zeros((l,20))
ii=0
dt=0.0001
N=120000000
Neq=20000000
Ndiff=N-Neq
repetitions=20
runs=500
T=runs*repetitions
for x in D:
	for y in range(1,21):
		times,avalues,jrvalues,j2values,state=[],[],[],[],[]
		file=open('/home/richard/outhome/time%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			avalues.append(float(row[0]))
			jrvalues.append(float(row[1]))
			j2values.append(float(row[2]))
			state.append(float(row[3]))
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
		if counts<9999:
			changespertime[ii][y-1]=counts
		else:
			changespertime[ii][y-1]=(counts/(avaluesa[9999]*repetitions+jrvaluesa[9999]+j2valuesa[9999]/Ndiff))*T		
#		plt.figure()		
#		plt.hist(intbt, bins=50)
#		plt.yscale('log')
#		plt.title("distribution of bursting time intervals")
#		plt.savefig('bdistajrj2%s%d%d.pdf' %(date,x,y))
#		plt.figure()
#		plt.hist(inteqt, bins=50)
#		plt.yscale('log')
#		plt.title("distribution of equilibrium time intervals")
#		plt.savefig('eqdistajrj2%s%d%d.pdf' %(date,x,y))
#		eqtot=np.sum(inteqt)
#		btot=np.sum(intbt)
#		eqrel=eqtot/(eqtot+btot)
#		brel=btot/(eqtot+btot)
#		eqtottime[ii][y-1]=eqtot/counteq
#		btottime[ii][y-1]=btot/countb
#		eqreltime[ii][y-1]=eqrel
#		breltime[ii][y-1]=brel
	ii=ii+1	


plt.figure()
xs2=np.arange(-0.175,0.325,0.025)
popt,pcov = curve_fit(deffana, xs2, changespertime[0,:])
params2[0]=popt[0]
params2[1]=popt[1]
params2[2]=popt[2]
params2[3]=popt[3]
params2[4]=popt[4]
params2[5]=popt[5]
eqfile2 = open('changesfit%s%d.txt' %(date,D[0]),'w')
for k4 in range(0,6): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.close() 
t=np.arange(-0.2,0.3,0.01)	
plt.xlabel('bias current')
plt.ylabel('changes')
plt.yscale('log')
for n in range(0,l):
	plt.plot(xs2,changespertime[n,:],label='D=%s' %D[n])
plt.plot(t,deffana(t,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),label='fit')
plt.legend()
plt.savefig('totalchanges%s%2.0f.pdf' %(date,D[0]))
