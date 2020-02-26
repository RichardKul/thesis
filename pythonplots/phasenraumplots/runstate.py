#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

v=open('/home/richard/cppworkshop/phasev11stime0.1.txt',"r")
vv=[]
for k in v:
 num = float(k.rstrip())
 vv.append(num)
av=np.array(vv)
vl=len(av)
n=open('/home/richard/cppworkshop/phasen11stime0.1.txt',"r")
vn=[]
for k2 in n:
 num = float(k2.rstrip())
 vn.append(num)
an=np.array(vn)
nl=len(an)
z=open('/home/richard/cppworkshop/phasez11stime0.1.txt',"r")
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

I=0.1
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
#plt.ylim(-0.1,an[-1])
plt.suptitle('equil. time [ms]')
plt.xlabel('initial $V_0$')
plt.ylabel('initial $n_0$')
#plt.legend()
plt.tight_layout()
plt.savefig('contoureqtime01wn3.pdf')
