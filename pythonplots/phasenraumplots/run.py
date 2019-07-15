#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

v=open('/home/richard/outhome/phasevneu4.txt',"r")
vv=[]
for k in v:
 num = float(k.rstrip())
 vv.append(num)
av=np.array(vv)
vl=len(av)
n=open('/home/richard/outhome/phasenneu4.txt',"r")
vn=[]
for k2 in n:
 num = float(k2.rstrip())
 vn.append(num)
an=np.array(vn)
nl=len(an)
z=open('/home/richard/outhome/phasezneu4.txt',"r")
vz=np.zeros((nl,vl))
v=0
for k3 in z:
 rowz=k3.split()
 rowza=np.array(rowz)
 for f in range(0,vl):
  vz[v][f]=rowza[f]
 v=v+1
t=np.arange(-70,-17,0.1)
vnc=(4-8*(t+80)-20*(t-60)/(1+np.exp((-20-t)/15)))/(9*(t+90))
nnc=1/(1+np.exp((-25-t)/5))
b=1.65+0.0261*t
#s=0.76*np.sin((t+56)*np.pi/80)

plt.figure()
plt.contourf(av,an,vz)
plt.colorbar()
plt.plot(t,vnc,label='v-nullcline')
plt.plot(t,nnc,label='n-nullcline')
plt.plot(t,b,label='barrier')
#plt.plot(t,s,label='barrier2')
plt.plot(-62.5701,0.00054509, 'bo')
plt.suptitle('spike count per starting point')
plt.xlabel('initial membrane voltage')
plt.ylabel('initial gating variable')
plt.savefig('contourneu4neuwn.pdf')
