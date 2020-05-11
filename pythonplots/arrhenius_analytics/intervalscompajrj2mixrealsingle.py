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
def fund(x,a,b,alpha):
	return a * (1/x)**alpha * np.exp(-b * x)

D1=[35]
D2=[30]
D3=[40,50]
Dvar=[30]
D=D1+D2+D3+Dvar
l1=len(D1)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l1+l2+l3+lvar
date2='new'+'realfast11jjem2sh'
date1='new'+'realfast19jjem2st'
date3='new'+'realfast11jjem2st'
datevar=['new'+'realfast11jjem2','new'+'realfast11jjem2sh','new'+'realfast11jjem2']
istart=9
ivalues=1
params2=np.zeros(6)
params=np.zeros((4,ivalues))
changespertime=np.zeros((l,ivalues))
eqtottime=np.zeros((l,ivalues))
eqreltime=np.zeros((l,ivalues))
btottime=np.zeros((l,ivalues))
breltime=np.zeros((l,ivalues))
intdelay=np.zeros((l,ivalues))

matplotlib.rcParams.update({'font.size': 22})

ii=0
dt=0.0005	
N=220000000
Neq=20000000
Ndiff=N-Neq
runs=500
repetitions=20
for x in D2:
	for y in range(istart,istart+ivalues):
		avalues,jrvalues,j2values,state=[],[],[],[]
		file=open('/home/richard/outhome/time%s%d%d.txt' % (date2,x,y),"r")
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
		intdelay[ii][y-istart]=sum(i < 500 for i in intbt)/countb
		counts=counteq+countb
		countsrel=counts/Ndiff
		changespertime[ii][y-istart]=countsrel/(repetitions*runs)
		plt.figure()	
		plt.xlabel('interval length [s]')
		plt.ylabel('number of intervals')	
		plt.hist(intbt/1000, bins=50)
		plt.yscale('log')
		plt.title("run. intervals, $I=%.2f$, $D=%.2f$" %((-5+y)*0.02,x/100), fontsize=22)
		plt.tight_layout()
		plt.savefig('bdistplotmaster2.pdf')# %(date1,x,y))
		plt.figure()
		plt.xlabel('interval length [s]')
		plt.ylabel('number of intervals')
		plt.hist(inteqt/1000, bins=50)
		plt.yscale('log')
		plt.title("eq. intervals, $I=%.2f$, $D=%.2f$" %((-5+y)*0.02,x/100), fontsize=22)
		plt.tight_layout()
		plt.savefig('eqdistplotmaster2.pdf')# %(date1,x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-istart]=eqtot/counteq
		btottime[ii][y-istart]=btot/countb
		eqreltime[ii][y-istart]=eqrel
		breltime[ii][y-istart]=brel
	ii=ii+1	
