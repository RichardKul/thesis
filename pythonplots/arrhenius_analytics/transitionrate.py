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
def funcrate(x,ap,am,bp,bm):
	return 2*am * np.exp(-bm * x)*ap * np.exp(-bp * x)/(am * np.exp(-bm * x)+ap * np.exp(-bp * x))

params2=np.zeros(6)
params=np.zeros((4,20))
eqtottime=np.zeros((3,20))
eqreltime=np.zeros((3,20))
btottime=np.zeros((3,20))
breltime=np.zeros((3,20))
ii=0
dt=0.00001

N2=100000000
Neq2=10000000
Ndiff2=N2-Neq2
runs=100
repetitions=100
for x in [30,40,50]:
	for y in range(1,21):
		times,avalues,jrvalues,j2values,state=[],[],[],[],[]
		file=open('/home/richard/outhome/time7m%d%d.txt' % (x,y),"r")
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
		eqtot=np.sum(inteqt)
		btot=np.sum(intbt)
		eqrel=eqtot/(eqtot+btot)
		brel=btot/(eqtot+btot)
		eqtottime[ii][y-1]=eqtot/counteq
		btottime[ii][y-1]=btot/countb
		eqreltime[ii][y-1]=eqrel
		breltime[ii][y-1]=brel
	ii=ii+1
	inteqt[c]=inteq[c]				

for k2 in range(0,20):
	plt.figure()
	xs=[1/3,1/4,1/5]
	xv=np.arange(0.1,0.5,0.01)
	popt,pcov = curve_fit(func, xs, 1/btottime[:,k2])
	params[0][k2]=popt[0]
	params[1][k2]=popt[1]
	ap=params[0][k2]
	bp=params[1][k2]
	popt,pcov = curve_fit(func, xs, 1/eqtottime[:,k2])
	params[2][k2]=popt[0]
	params[3][k2]=popt[1]
	am=params[2][k2]
	bm=params[3][k2]
	eqfile = open('paramsh16m%d.txt' % (k2),'w')
	for k3 in range(0,4): 
		eqfile.write('%.6f\n'%params[k3][k2]) 
	eqfile.close() 
	plt.figure()
	plt.xlabel('inverse noise intensity 1/D')
	plt.ylabel('total transition rate')
	plt.yscale('log')
	plt.plot(xs,(2*(1/btottime[:,k2])*(1/eqtottime[:,k2]))/(1/eqtottime[:,k2]+1/btottime[:,k2]),'bo',label='simulation')
	plt.plot(xv,funcrate(xv,ap,am,bp,bm),label='theory')
	plt.savefig('cortimeshort16m%d.pdf' %(k2))
plt.figure()
xold=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current I')
plt.ylabel('correlation time')
plt.yscale('log')
plt.plot(xold,1/(1/btottime[0,:]+1/eqtottime[0,:]),label='D=3')
plt.plot(xold,1/(1/btottime[1,:]+1/eqtottime[1,:]),label='D=4')
plt.plot(xold,1/(1/btottime[2,:]+1/eqtottime[2,:]),label='D=5')
plt.savefig('altcortimeshort16m.pdf')

plt.figure()
xnew=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current')
plt.ylabel('prefactor')
plt.plot(xnew,params[0,:],label='burst to eq')
plt.plot(xnew,params[2,:],label='eq to burst')
plt.legend()
plt.savefig('prefacshort16m.pdf')
plt.figure()
plt.xlabel('bias current')
plt.ylabel('potential barrier')
plt.plot(xnew,params[1,:],label='burst to eq')
plt.plot(xnew,params[3,:],label='eq to burst')
plt.plot(xnew,2*params[1,:],label='2x burst to eq')
plt.plot(xnew,2*params[3,:],label='2x eq to burst')
plt.legend()
plt.savefig('barriershort216m.pdf')
eqfile3 = open('barrierex16m.txt','w')
for k3 in range(0,11): 
	eqfile3.write('%.6f %6.f %6.f %6.f %6.f\n'%(xnew[k3],params[0][k3],params[1][k3],params[2][k3],params[3][k3])) 
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
eqfile2 = open('parambarriershort16m.txt','w')
for k4 in range(0,6): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.close() 

t=np.arange(0.5,3,0.01)
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
plt.savefig('barriercompshort16m.pdf')


