#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def comp(x,a,b):
	return a*x+b

file=open('/home/richard/mastergit/NetBeansProjects/detmodelreal/countanhopfa16.txt',"r")
col,colx=[],[]
for k in file:
	row=k.split()
	colx.append(float(row[0]))
	col.append(float(row[1]))
cola=np.array(col)
colxa=np.array(colx)
a=174/9.8
a2=82/4

t=np.arange(0,10,0.1)
xs=np.arange(0,10,0.2)
plt.xlabel('timestep')
plt.ylabel('spike count')
#plt.xlim(0,4)
#plt.ylim(1.2,1.4)
plt.yscale('log')
plt.xscale('log')
plt.plot(colxa,cola)
#plt.plot(t,comp(t,a2,602)/500,label='linear appr')
#plt.legend()
plt.tight_layout()
plt.savefig('detmotimeanhopf4.pdf')


