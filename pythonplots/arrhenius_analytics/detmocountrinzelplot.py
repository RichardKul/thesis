#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def comp(x,a,b):
	return a/x+b

D=np.array([25,30,40,50])

params=4
dvalues=6
b=np.zeros((params,dvalues))
ii=0
eqfile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/detmocountparam5.txt','r')
for k in eqfile:
	col=[]
	row=k.split()
	for lll in range(0,params):	
		col.append(float(row[lll]))
	cola=np.array(col)
	for kk in range(0,params):
		b[kk][ii]=cola[kk]
	ii=ii+1
eqfile.close() 
a1=np.zeros(params)
a2=np.zeros(params)
eqfile = open('detmocountparamparam2.txt','w')
for n in range(0,params):
	popt,pcov = curve_fit(comp, D, b[n,2:6])
	a1[n]=popt[0]
	a2[n]=popt[1]
	for k3 in range(0,2): 
		eqfile.write('%.6f '%popt[k3])
	eqfile.write('\n') 
eqfile.close() 

t=np.arange(20,50,0.1)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.xlim(0,4)
#plt.ylim(1.2,1.4)
#plt.yscale('log')
for n in range(0,params):
	plt.figure()
	plt.plot(D,b[n,2:6])
	plt.plot(t,comp(t,a1[n],a2[n]))
	plt.savefig('detmocountrinzelparam2%.0f.pdf' %n)
for n in range(0,params):
	plt.figure()
	plt.plot(1/t,comp(t,a1[n],a2[n]))
	plt.plot(1/D,b[n,2:6])
	plt.savefig('detmocountrinzelparaminv2%.0f.pdf' %n)


