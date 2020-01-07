#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

date3='raterealrinzel25orealrinzel15ninv0'
date2=''

istart=1
ivalues=11
l=3

D=[200,300,500]
Da=np.array(D)
btoeq=np.zeros((l,ivalues))
eqtob=np.zeros((l,ivalues))

for k2 in range(0,ivalues):
	x=[]
	ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/%s%d.txt' %(date3+date2,k2),'r')
	for k4 in ratefile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	for k in range(0,l):
		btoeq[k][k2]=1/ax[k]
		eqtob[k][k2]=1/ax[k+l]

xnew=np.arange(-22.5+istart,-22.5+istart+ivalues)*0.8
#xnew=np.arange(-5+istart,-5+istart+ivalues)*0.02
colorv=['y','g','b','r','c','k']
plt.figure()
plt.yscale('log')
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('transition rate w $[10^3s^{-1}]$')
for n in range(0,l):
	plt.plot(xnew,btoeq[n,:],colorv[n]+'x',label='D=%.0f, burst to eq'%(Da[n]*0.1))
	plt.plot(xnew,btoeq[n,:],colorv[n])
	plt.plot(xnew,eqtob[n,:],colorv[n]+'+')
	plt.plot(xnew,eqtob[n,:],colorv[n])
plt.legend()
plt.savefig('tranrates.pdf')
