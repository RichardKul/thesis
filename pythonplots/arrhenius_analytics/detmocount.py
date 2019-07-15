#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def comp(x,a,b):
	return a*x+b

file=open('/home/richard/NetBeansProjects/detmodel/countI.txt',"r")
col=[]
for k in file:
	row=k.split()
	col.append(float(row[1]))
cola=np.array(col)
a=174/9.8
a2=82/4

t=np.arange(0,10,0.1)
xs=np.arange(0,10,0.2)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
plt.xlim(0,4)
plt.ylim(1.2,1.4)
#plt.yscale('log')
plt.plot(xs,cola/500,label='measured')
plt.plot(t,comp(t,a2,602)/500,label='linear appr')
plt.legend()
plt.savefig('detmocount04.pdf')


