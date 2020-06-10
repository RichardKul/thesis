#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def comp(x,b,c,d,e):
	return b*x**3+c*x**2+d*x+e

matplotlib.rcParams.update({'font.size': 20})

N=5000000
dt=0.005
T=N*dt
dvalues=6
ivalues=51
D=[150,200,250,300,400,500]
count=np.zeros((dvalues,ivalues))
ii=0
file=open('/home/richard/NetBeansProjects/detrinzel/countrinzelnoise5.txt',"r")
colx=[]
for k in file:
	col=[]
	row=k.split()
	colx.append(float(row[0]))
	for l in range(1,dvalues+1):	
		col.append(float(row[l]))
	cola=np.array(col)
	for l in range(0,dvalues):
		count[l][ii]=cola[l]
	ii=ii+1
colxa=np.array(colx)
a=174/9.8
a2=82/4
#a=np.zeros(dvalues)
b=np.zeros(dvalues)
c=np.zeros(dvalues)
d=np.zeros(dvalues)
e=np.zeros(dvalues)
#xnew=np.arange(-0.19,0.31,0.01)
eqfile = open('detmocountparam6.txt','w')
for n in range(0,dvalues):
	popt,pcov = curve_fit(comp, colxa, count[n,:]/T)
#	a[n]=popt[0]
	b[n]=popt[0]
	c[n]=popt[1]
	d[n]=popt[2]
	e[n]=popt[3]
	for k3 in range(0,4): 
		eqfile.write('%.6f '%popt[k3])
	eqfile.write('\n') 
eqfile.close() 
#colorv=['y','g','b','r','c','k']
colorv=[ '#1f77b4', '#ff7f0e', '#2ca02c','#d62728','#9467bd']
t=np.arange(-21,-6,0.1)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('spiking firing rate $v_0$ [$s^{-1}$]')
#plt.xlim(0,4)
#plt.ylim(1.2,1.4)
#plt.yscale('log')
for n in range(2,dvalues):
	plt.plot(colxa,count[n,:]/T*1000,'x',color=colorv[n-1],label='D=%.2s' %(D[n]/10))
	plt.plot(t,comp(t,b[n],c[n],d[n],e[n])*1000,colorv[n-1])
#plt.plot(t,comp(t,popt[0],popt[1]),label='linear appr %f %f'%(popt[0],popt[1]))
plt.legend()
plt.tight_layout()
plt.savefig('detmocountrinzelcomp6big.pdf')


