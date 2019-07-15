#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def deff(r0p,r0m,ap,am,bp,bm,t,D,v):
	return (v**2*r(r0p,ap,bp,t,D)*r(r0m,am,bm,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**3)

file=open('/home/richard/NetBeansProjects/parambarrier7m.txt',"r")
x=[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
ax=np.array(x)

t=np.arange(0,4,0.01)
ap=ax[0]
am=ax[3]
bp=ax[1]
bm=ax[4]
r0p=ax[2]
r0m=ax[5]
v=1.33

plt.figure()
plt.xlabel('bias current I')
plt.ylabel('barrier height')
#plt.yscale('log')
plt.plot(t,barrier(am,bm,t),label='eq to burst')
plt.plot(t,barrier(ap,bp,t),label='burst to eq')
plt.plot(t,2*barrier(am,bm,t),label='2x eq to burst')
plt.plot(t,2*barrier(ap,bp,t),label='2x burst to eq')
plt.legend()
plt.savefig('calcbarr7m.pdf')
