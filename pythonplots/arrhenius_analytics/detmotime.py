#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 20})

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
fig, ax=plt.subplots()
plt.xlabel('timestep [ms]')
plt.ylabel('spike count')
#plt.xlim(0,4)
#plt.ylim(1.2,1.4)
#plt.yscale('log')
plt.xscale('log')
plt.plot(colxa,cola)
#plt.plot(t,comp(t,a2,602)/500,label='linear appr')
#plt.legend()
ax.yaxis.set_major_locator(plt.MaxNLocator(3))
plt.tight_layout()
plt.savefig('detmotimeanhopfbig2.pdf')


