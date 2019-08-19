#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def comp(x,a,b):
	return a*x+b

N=50000000
dt=0.0005
T=N*dt
file=open('/home/richard/mastergit/NetBeansProjects/detmodel/countI9a.txt',"r")
col=[]
for k in file:
	row=k.split()
	col.append(float(row[1]))
cola=np.array(col)
a=174/9.8
a2=82/4

xnew=np.arange(-0.19,0.31,0.01)
popt,pcov = curve_fit(comp, xnew, cola/T)

t=np.arange(-0.2,0.3,0.01)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.xlim(0,4)
#plt.ylim(1.2,1.4)
#plt.yscale('log')
plt.plot(xnew,cola/T,label='measured')
plt.plot(t,comp(t,popt[0],popt[1]),label='linear appr %f %f'%(popt[0],popt[1]))
plt.legend()
plt.savefig('detmocount9a.pdf')


