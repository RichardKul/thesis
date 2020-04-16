#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

date='anhopfa2sh'
I=46

v=open('/home/richard/cppworkshop/phasev%s%d.txt' %(date,I),"r")
vv=[]
for k in v:
 num = float(k.rstrip())
 vv.append(num)
av=np.array(vv)
vl=len(av)
n=open('/home/richard/cppworkshop/phasen%s%d.txt' %(date,I),"r")
vn=[]
for k2 in n:
 num = float(k2.rstrip())
 vn.append(num)
an=np.array(vn)
nl=len(an)
z=open('/home/richard/cppworkshop/phasez%s%d.txt' %(date,I),"r")
vz=np.zeros((nl,vl))
v=0
for k3 in z:
 rowz=k3.split()
 rowza=np.array(rowz)
 for f in range(0,vl):
  vz[v][f]=rowza[f]
 v=v+1

def finf(x,v12,k):
	return 1/(1+np.exp((v12-x)/k))
def vnc(x,I,v12,k,gL,EL,gNa,ENa,gK,EK):
	return (I-gL*(x-EL)-gNa*finf(x,v12,k)*(x-ENa))/(gK*(x-EK)) 


gL=1
EL=-78
gNa=4
ENa=60
gK=4
EK=-90
km=7 
vm=-30
kn=5 
vn=-45
tau=1 
#b=1.65+0.0261*t
#s=0.76*np.sin((t+56)*np.pi/80)
t=np.arange(av[0],av[-1]+0.01,0.01)

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
plt.contourf(av,an,abs(vz))
plt.colorbar()
plt.plot(t,vnc(t,I,vm,km,gL,EL,gNa,ENa,gK,EK),label='v-nullcline')
plt.plot(t,finf(t,vn,kn),label='n-nullcline')
#plt.plot(t,b,label='barrier')
#plt.plot(t,s,label='barrier2')
#plt.plot(-62.5701,0.00054509, 'bo')
plt.ylim(-0.1,an[-1])
#plt.suptitle('equil. time [ms]')
plt.suptitle('spikes/starting point')
plt.xlabel('initial $V_0$')
plt.ylabel('initial $n_0$')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('contour%s.pdf' %date)
