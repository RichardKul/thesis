#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


def an(V):
	return (10-V)/(100*(np.exp((10-V)/10)-1))
def ah(V):
	return 7*np.exp(-V/20)/100
def bn(V):
	return np.exp(-V/80)/8
def bh(V):
	return 1/(np.exp((30-V)/10)+1)
def hinf(V):
	return ah(V)/(ah(V)+bh(V))
def ninf(V):
	return an(V)/(an(V)+bn(V))
def winf(S,V):
	return S*(ninf(V)+S*(1-hinf(V)))/(1+S**2) # entspricht w-isokline

gL=0.3
VL=10
gNa=120
VNa=115
gK=36
VK=12 #12
phi=3.82
Iapp=-10 #-12
SV=1.27

vec=np.arange(-50,150,0.01)

vrinzelfile = open('vncrinzel%d.txt' %Iapp, "r")
colx,col=[],[]
for k in vrinzelfile:
	row=k.split()
	colx.append(float(row[0]))
	col.append(float(row[1]))
edgepointfile= open('edgerinzel%d.txt' %Iapp, "r")
edge=[]
for k2 in edgepointfile:
	row2=k2.split()
	edge.append(int(row2[0]))


date='rinzela12time'
I=-10

v=open('/home/richard/cppworkshop/phasev%s%d.txt' %(date,I),"r")
#v=open('/home/richard/outhome/phasev%s%d.txt' %(date,I),"r")
vv=[]
for k in v:
 num = float(k.rstrip())
 vv.append(num)
av=np.array(vv)
vl=len(av)
n=open('/home/richard/cppworkshop/phasen%s%d.txt' %(date,I),"r")
#n=open('/home/richard/outhome/phasen%s%d.txt' %(date,I),"r")
vn=[]
for k2 in n:
 num = float(k2.rstrip())
 vn.append(num)
avn=np.array(vn)
nl=len(avn)
z=open('/home/richard/cppworkshop/phasez%s%d.txt' %(date,I),"r")
#z=open('/home/richard/outhome/phasez%s%d.txt' %(date,I),"r")

vz=np.zeros((nl,vl))
v=0
for k3 in z:
 rowz=k3.split()
 rowza=np.array(rowz)
 for f in range(0,vl):
  vz[v][f]=rowza[f]
 v=v+1

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
plt.contourf(av,avn,abs(vz))
plt.colorbar()
plt.plot(colx[edge[0]:edge[1]],col[edge[0]:edge[1]],'blue',label='V-nullcline')
plt.plot(colx[edge[2]:edge[3]],col[edge[2]:edge[3]],'blue')
plt.plot(vec,winf(SV,vec),'orange',label='W-nullcline')
#plt.plot(t,b,label='barrier')
#plt.plot(t,s,label='barrier2')
#plt.plot(-62.5701,0.00054509, 'bo')
plt.ylim(-0.1,avn[-1])
#plt.xlim(8,12)
plt.suptitle('equil. time [ms]')
#plt.suptitle('spikes/starting point')
plt.xlabel('initial $V_0$')
plt.ylabel('initial $W_0$')
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig('contour%s.pdf' %date)
