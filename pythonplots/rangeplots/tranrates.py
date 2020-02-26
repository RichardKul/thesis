#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

date3='ratenewrealfast11jjem2shnewrealfast19jjem2st'
date2=''

istart=1
ivalues=20
l=5

D=[35,45,40,50,30]
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
xnew=np.arange(-5+istart,-5+istart+ivalues)*0.02
#colorv=['y','g','b','r','c','k']
colorv=['b','r','k','g','y']
plt.figure()
plt.yscale('log')
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('transition rate r $[s^{-1}]$')
for n in range(0,2):
	plt.plot(xnew,btoeq[n,:],colorv[n]+'x',label='D=%.2f, run to eq'%(Da[n]*0.01))
	plt.plot(xnew,btoeq[n,:],colorv[n])
	plt.plot(xnew,eqtob[n,:],colorv[n]+'+')
	plt.plot(xnew,eqtob[n,:],colorv[n])
for n in range(4,5):
	plt.plot(xnew,btoeq[n,:],colorv[n]+'x',label='D=%.2f, run to eq'%(Da[n]*0.01))
	plt.plot(xnew,btoeq[n,:],colorv[n])
	plt.plot(xnew,eqtob[n,:],colorv[n]+'+')
	plt.plot(xnew,eqtob[n,:],colorv[n])
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('tranratesneur3.pdf')
