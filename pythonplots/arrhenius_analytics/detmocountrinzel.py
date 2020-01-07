#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def comp(x,b,c,d,e):
	return b*x**3+c*x**2+d*x+e

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
eqfile = open('detmocountparam5.txt','w')
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
colorv=['y','g','b','r','c','k']
t=np.arange(-21,-6,0.1)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.xlim(0,4)
#plt.ylim(1.2,1.4)
#plt.yscale('log')
for n in range(0,dvalues):
	plt.plot(colxa,count[n,:]/T,colorv[n]+'o',label='D=%s' %(D[n]/10))
	plt.plot(t,comp(t,b[n],c[n],d[n],e[n]),colorv[n])
#plt.plot(t,comp(t,popt[0],popt[1]),label='linear appr %f %f'%(popt[0],popt[1]))
plt.legend()
plt.savefig('detmocountrinzelcomp5.pdf')


