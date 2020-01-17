#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def an(V):
	return (10-V)/(100*(np.exp((10-V)/10)-1))
def am(V):
	return (25-V)/(10*(np.exp((25-V)/10)-1))
def ah(V):
	return 7*np.exp(-V/20)/100
def bn(V):
	return np.exp(-V/80)/8
def bm(V):
	return 4*np.exp(-V/18)
def bh(V):
	return 1/(np.exp((30-V)/10)+1)
def minf(V):
	return am(V)/(am(V)+bm(V))
def hinf(V):
	return ah(V)/(ah(V)+bh(V))
def ninf(V):
	return an(V)/(an(V)+bn(V))
def tau(V):
	return 5*np.exp(-(((V+10)/55)**2))+1
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

def viso(W,V):
	return W**4+(Iapp-gNa*((minf(V))**3)*(1-W)*(V-VNa)-gL*(V-VL))/(-gK*((1/SV)**4)*(V-VK))

v=np.arange(-50,150,0.01)
lv=len(v)
w=np.arange(0,1,0.0001)

vvalues=np.zeros(lv)

eqpoifile = open('eqpointsrinzelcorI%d.txt' %Iapp, 'w')

for k in range(0,lv):
	vnew=w[np.where(abs(viso(w,v[k]))==min(abs(viso(w,v[k]))))[0][0]]
	if abs(viso(vnew,v[k]))<0.01:
		vvalues[k]=vnew
	else:
		vvalues[k]=0
	if k>0:
		if (winf(SV,v[k])-vold)*(winf(SV,v[k])-vnew)<0:
			eqpoifile.write('%.2f\n'%v[k]) 	
	vold=vnew
eqpoifile.close() 
plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('recovery variable W')
plt.plot(v,winf(SV,v),label='W-nullcline')
plt.plot(v,vvalues,label='V-nullcline')
plt.legend()
plt.savefig('rinzelclinescor10.pdf')

