#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def finf(x,v12,k):
	return 1/(1+np.exp((v12-x)/k))
def vnc(x,I,v12,k,gL,EL,gNa,ENa,gK,EK):
	return (I-gL*(x-EL)-gNa*finf(x,v12,k)*(x-ENa))/(gK*(x-EK)) 

I=0.2
gL=0.3
EL=-80
gNa=1
ENa=60
gK=0.4
EK=-90
km=14 
vm=-18
kn=5 
vn=-25
tau=3 

file=open('/home/richard/mastergit/NetBeansProjects/inapik/inapreali20.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)

t=np.arange(-55,5,0.1)

plt.figure()
plt.xlabel('membrane voltage V [mV]')
plt.ylabel('gating variable n')
plt.plot(t,vnc(t,I,vm,km,gL,EL,gNa,ENa,gK,EK),label='v-nullcline')
plt.plot(t,finf(t,vn,kn),label='n-nullcline')
#plt.gcf().subplots_adjust(left=0.15)
plt.plot(ax,ay)
plt.legend()
plt.savefig('inapreali20wn.pdf')
