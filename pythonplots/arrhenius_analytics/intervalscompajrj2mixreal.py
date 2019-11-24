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

D1=[35]
D2=[45]
D3=[40,50]
Dvar=[30]
D=D1+D2+D3+Dvar
l1=len(D1)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l1+l2+l3+lvar
date1='new'+'realfast11jjem2sh'
date2='new'+'realfast19jjem2st'
date3='new'+'realfast11jjem2st'
datevar=['new'+'realfast11jjem2','new'+'realfast11jjem2sh','new'+'realfast11jjem2']
istart=1
ivalues=20
params2=np.zeros(6)
params=np.zeros((4,ivalues))
changespertime=np.zeros((l,ivalues))
eqtottime=np.zeros((l,ivalues))
eqreltime=np.zeros((l,ivalues))
btottime=np.zeros((l,ivalues))
breltime=np.zeros((l,ivalues))
ii=0
dt=0.0005	
N=220000000
Neq=20000000
Ndiff=N-Neq
runs=500
repetitions=20
for x in D1:
	for y in range(istart,istart+ivalues):
		avalues,jrvalues,j2values,state=[],[],[],[]
		file=open('/home/richard/outhome/time%s%d%d.txt' % (date1,x,y),"r")
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
		changespertime[ii][y-istart]=countsrel/(repetitions*runs)
		plt.figure()	
		plt.xlabel('interval length [ms]')
		plt.ylabel('number of intervals')	
		plt.hist(intbt, bins=50)
		plt.yscale('log')
		plt.title("distribution of bursting time intervals")
		plt.savefig('bdistajrj2%s%d%d.pdf' %(date1,x,y))
		plt.figure()
		plt.xlabel('interval length [ms]')
		plt.ylabel('number of intervals')
		plt.hist(inteqt, bins=50)
		plt.yscale('log')
		plt.title("distribution of equilibrium time intervals")
		plt.savefig('eqdistajrj2%s%d%d.pdf' %(date1,x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-istart]=eqtot/counteq
		btottime[ii][y-istart]=btot/countb
		eqreltime[ii][y-istart]=eqrel
		breltime[ii][y-istart]=brel
	ii=ii+1	

	
N=120000000
Ndiff=N-Neq

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
		
		counts=counteq+countb
		countsrel=counts/Ndiff
		changespertime[ii][y-istart]=countsrel/(repetitions*runs)
		plt.figure()		
		plt.xlabel('interval length [ms]')
		plt.ylabel('number of intervals')	
		plt.hist(intbt, bins=50)
		plt.yscale('log')
		plt.title("distribution of bursting time intervals")
		plt.savefig('bdistajrj2%s%d%d.pdf' %(date2,x,y))
		plt.figure()
		plt.xlabel('interval length [ms]')
		plt.ylabel('number of intervals')	
		plt.hist(inteqt, bins=50)
		plt.yscale('log')
		plt.title("distribution of equilibrium time intervals")
		plt.savefig('eqdistajrj2%s%d%d.pdf' %(date2,x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-istart]=eqtot/counteq
		btottime[ii][y-istart]=btot/countb
		eqreltime[ii][y-istart]=eqrel
		breltime[ii][y-istart]=brel
	ii=ii+1	


for x in D3:
	for y in range(istart,istart+ivalues):
		avalues,jrvalues,j2values,state=[],[],[],[]
		file=open('/home/richard/outhome/time%s%d%d.txt' % (date3,x,y),"r")
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
		changespertime[ii][y-istart]=countsrel/(repetitions*runs)
		plt.figure()		
		plt.xlabel('interval length [ms]')
		plt.ylabel('number of intervals')	
		plt.hist(intbt, bins=50)
		plt.yscale('log')
		plt.title("distribution of bursting time intervals")
		plt.savefig('bdistajrj2%s%d%d.pdf' %(date3,x,y))
		plt.figure()
		plt.xlabel('interval length [ms]')
		plt.ylabel('number of intervals')	
		plt.hist(inteqt, bins=50)
		plt.yscale('log')
		plt.title("distribution of equilibrium time intervals")
		plt.savefig('eqdistajrj2%s%d%d.pdf' %(date3,x,y))
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-istart]=eqtot/counteq
		btottime[ii][y-istart]=btot/countb
		eqreltime[ii][y-istart]=eqrel
		breltime[ii][y-istart]=brel
	ii=ii+1	

dt=0.0005	
repetitions=20
runs=500
yvar=[4,13,3]
yvalues=len(yvar)
Nvar=[520000000,220000000,520000000]
Neq=20000000
for x in Dvar:
	for kk in range(0,yvalues):
		Ndiff=Nvar[kk]-Neq
		ystart=1
		jj=0
		while jj < kk:
			ystart=ystart+yvar[jj]
			jj=jj+1
		for y in range(ystart,ystart+yvar[kk]):
			avalues,jrvalues,j2values,state=[],[],[],[]
			file=open('/home/richard/outhome/time%s%d%d.txt' % (datevar[kk],x,y),"r")
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
					if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/	Ndiff>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff:
						inteq[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-jrvaluesa[2*z])*Ndiff+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
						counteq=counteq+1
					else:
						break	
				inteqt=np.zeros(counteq)
				for c in range(0,counteq):
					inteqt[c]=inteq[c]				
				for z in range(0,4999):
					if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/	Ndiff>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff:
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
					if avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/	Ndiff>avaluesa[2*z]*repetitions+jrvaluesa[2*z]+j2valuesa[2*z]/Ndiff:
						intb[z]=(((avaluesa[2*z+1]-avaluesa[2*z])*repetitions+jrvaluesa[2*z+1]-	jrvaluesa[2*z])*Ndiff+j2valuesa[2*z+1]-j2valuesa[2*z])*dt
						countb=countb+1
					else:
						break	
				intbt=np.zeros(countb)
				for c in range(0,countb):
					intbt[c]=intb[c]				
				for z in range(0,4999):
					if avaluesa[2*z+2]*repetitions+jrvaluesa[2*z+2]+j2valuesa[2*z+2]/	Ndiff>avaluesa[2*z+1]*repetitions+jrvaluesa[2*z+1]+j2valuesa[2*z+1]/Ndiff:
						inteq[z]=(((avaluesa[2*z+2]-avaluesa[2*z+1])*repetitions+jrvaluesa[2*z+2]-jrvaluesa[2*z+1])*Ndiff+j2valuesa[2*z+2]-j2valuesa[2*z+1])*dt
						counteq=counteq+1
					else:
						break
				inteqt=np.zeros(counteq)
				for c in range(0,counteq):
					inteqt[c]=inteq[c]	
			
			counts=counteq+countb
			countsrel=counts/Ndiff
			changespertime[ii][y-istart]=countsrel/(repetitions*runs)
			plt.figure()		
			plt.xlabel('interval length [ms]')
			plt.ylabel('number of intervals')	
			plt.hist(intbt, bins=50)
			plt.yscale('log')
			plt.title("distribution of bursting time intervals")
			plt.savefig('bdistajrj2%s%d%d.pdf' %(datevar[0],x,y))
			plt.figure()
			plt.xlabel('interval length [ms]')
			plt.ylabel('number of intervals')	
			plt.hist(inteqt, bins=50)
			plt.yscale('log')
			plt.title("distribution of equilibrium time intervals")
			plt.savefig('eqdistajrj2%s%d%d.pdf' %(datevar[0],x,y))
			eqtot=np.sum(inteqt)
			btot=np.sum(intbt)
			eqrel=eqtot/(eqtot+btot)
			brel=btot/(eqtot+btot)
			eqtottime[ii][y-istart]=eqtot/counteq
			btottime[ii][y-istart]=btot/countb
			eqreltime[ii][y-istart]=eqrel
			breltime[ii][y-istart]=brel
	ii=ii+1	
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
if l > 1:
	for k2 in range(0,ivalues):
		plt.figure()
		xs=np.zeros(l)
		for xf in range(0,l):
			xs[xf]=100/D[xf]
		plt.suptitle('I=%.2f$\mu A/cm^2$' %(-0.1+0.02*k2))
		plt.xlabel('inverse noise intensity 1/D')
		plt.ylabel('transition rate w $[10^3s^{-1}]$')
		plt.yscale('log')
		plt.plot(xs,1/btottime[:,k2],'bo',label='burst to eq')
#plt.plot(xs,breltime[1,:],label='D=3,burst')
#plt.plot(xs,breltime[2,:],label='D=4,burst')
#plt.plot(xs,eqreltime[3,:],label='D=4,burst')
		plt.plot(xs,1/eqtottime[:,k2],'ro',label='eq to burst')
#plt.plot(xs,eqreltime[1,:],label='D=3,equilibrium')
#plt.plot(xs,eqreltime[2,:],label='D=4,equilibrium')
#plt.plot(xs,breltime[3,:],label='D=4,equilibrium')
		popt,pcov = curve_fit(func, xs, 1/btottime[:,k2])
		plt.plot(np.array(xs), func(np.array(xs), *popt), 'b-',label='fit burst to eq: r_0=%5.3f, U_+=%5.3f' % tuple(popt))
		params[0][k2]=popt[0]
		params[1][k2]=popt[1]
		popt,pcov = curve_fit(func, xs, 1/eqtottime[:,k2])
		plt.plot(np.array(xs), func(np.array(xs), *popt), 'r-',label='fiteq to burst: r_0=%5.3f, U_-=%5.3f' % tuple(popt))
		params[2][k2]=popt[0]
		params[3][k2]=popt[1]
		plt.legend()
		plt.savefig('arrheniustot%sfit%d.pdf' %(date1+date2,k2))
		eqfile = open('param%s%d.txt' % (date1+date2,k2),'w')
		for k3 in range(0,4): 
			eqfile.write('%.6f\n'%params[k3][k2]) 
		eqfile.close() 
		ratefile = open('rate%s%d.txt' %(date1+date2,k2),'w')
		for k4 in range(0,l):
			ratefile.write('%.6f\n'%btottime[k4][k2])
		for k4 in range(0,l):  
			ratefile.write('%.6f\n'%eqtottime[k4][k2]) 
		ratefile.close() 
		plt.figure()
		plt.xlabel('inverse noise intensity 1/D')
		plt.ylabel('correlation time')
		plt.yscale('log')
		plt.plot(xs,1/(1/btottime[:,k2]+1/eqtottime[:,k2]))
		plt.savefig('cortime%s%d.pdf' %(date1+date2,k2))
plt.figure()
xold=np.arange(-5+istart,-5+istart+ivalues)*0.02
plt.xlabel('bias current I')
plt.ylabel('correlation time')
plt.yscale('log')
for n in range(0,l1):
	plt.plot(xold,1/(1/btottime[n,:]+1/eqtottime[n,:]),label='D=%f' %(D1[n]*0.01))
for n in range(l1,l1+l2):
	plt.plot(xold,1/(1/btottime[n,:]+1/eqtottime[n,:]),label='D=%f' %(D2[n-l1]*0.01))
plt.savefig('altcortime%s.pdf' %(date1+date2))

plt.figure()
xnew=np.arange(-5+istart,-5+istart+ivalues)*0.02
plt.xlabel('bias current')
plt.ylabel('prefactor')
plt.plot(xnew,params[0,:],label='burst to eq')
plt.plot(xnew,params[2,:],label='eq to burst')
plt.legend()
plt.savefig('prefac%s.pdf' %(date1+date2))
plt.figure()
plt.xlabel('bias current')
plt.ylabel('potential barrier')
plt.plot(xnew,params[1,:],label='burst to eq')
plt.plot(xnew,params[3,:],label='eq to burst')
plt.plot(xnew,2*params[1,:],label='2x burst to eq')
plt.plot(xnew,2*params[3,:],label='2x eq to burst')
plt.legend()
plt.savefig('barrier2%s.pdf' %(date1+date2))
eqfile3 = open('barrierex%s.txt' %(date1+date2),'w')
for k3 in range(0,ivalues): 
	eqfile3.write('%.6f %.6f %.6f %.6f %.6f\n'%(xnew[k3],params[0][k3],params[1][k3],params[2][k3],params[3][k3])) 
eqfile3.close() 

popt,pcov = curve_fit(func2, xnew, params[1,:])
params2[0]=popt[0]
params2[1]=popt[1]
popt,pcov = curve_fit(func2, xnew, params[3,:])
params2[3]=popt[0]
params2[4]=popt[1]
popt,pcov = curve_fit(func3, xnew, params[0,:])
params2[2]=popt[0]
popt,pcov = curve_fit(func3, xnew, params[2,:])
params2[5]=popt[0]
eqfile2 = open('parambarrier%s.txt' %(date1+date2),'w')
for k4 in range(0,6): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.close() 

t=np.arange(-0.1,0.3,0.001)
plt.figure()
plt.xlabel('bias current')
plt.ylabel('potential barrier')
plt.plot(t,func2(t,params2[0],params2[1]),'y')
plt.plot(t,func2(t,params2[3],params2[4]),'y')
plt.plot(t,2*func2(t,params2[0],params2[1]),'y')
plt.plot(t,2*func2(t,params2[3],params2[4]),'y')
plt.plot(xnew,params[1,:],label='burst to eq')
plt.plot(xnew,params[3,:],label='eq to burst')
plt.plot(xnew,2*params[1,:],label='2x burst to eq')
plt.plot(xnew,2*params[3,:],label='2x eq to burst')

plt.legend()
plt.savefig('barriercomp%s.pdf' %(date1+date2))

#plt.figure()
#xs1=np.arange(-0.75,4.25,0.25)
#xs2=np.arange(0.2,4.2,0.2)	
#plt.xlabel('bias current')
#plt.ylabel('changes over time steps')
#plt.yscale('log')
#plt.plot(xs2,changespertime[0,:],label='D=1.2')
#plt.plot(xs1,changespertime[1,:],label='D=2')
#plt.plot(xs2,changespertime[2,:],label='D=3')
#plt.plot(xs2,changespertime[3,:],label='D=4')
#plt.plot(xs2,changespertime[4,:],label='D=5')
#plt.legend()
#plt.savefig('changespertime.pdf')
