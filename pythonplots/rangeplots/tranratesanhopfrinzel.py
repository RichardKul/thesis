#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

date3='raterealrinzelrange26d1'
date2='realrinzelrange26d1'



istart=1
ivalues=10
l=4

D=[250,300,400,500]
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
		btoeq[k][k2]=1000/ax[k]
		eqtob[k][k2]=1000/ax[k+l]

#xnew=np.arange(-22.5+istart,-22.5+istart+ivalues)*0.8
xnew=np.arange(-21.25+istart,-21.25+istart+ivalues)*0.8
#colorv=['y','g','b','r','c','k']
colorv=[ '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
#colorv=['b','r','k','g','y']
plt.figure()
plt.yscale('log')
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('transition rate r $[s^{-1}]$')
for n in range(0,l):
	plt.plot(xnew,btoeq[n,:],colorv[n+1],marker='x',label='D=%.0f, run to eq'%(Da[n]*0.1))
	plt.plot(xnew,btoeq[n,:],colorv[n+1])
	plt.plot(xnew,eqtob[n,:],colorv[n+1],marker='+')
	plt.plot(xnew,eqtob[n,:],colorv[n+1])
#plt.legend()
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('tranratesrinzel.pdf')
