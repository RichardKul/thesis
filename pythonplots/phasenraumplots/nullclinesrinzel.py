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
Iapp=-15 #-12
SV=1.27

def viso(W,V):
	return W**4+(Iapp-gNa*((minf(V))**3)*(1-W)*(V-VNa)-gL*(V-VL))/(-gK*((1/SV)**4)*(V-VK))

v=np.arange(-50,150,0.01)
lv=len(v)
w=np.arange(-0.1,1.1,0.0001)

vvalues=np.zeros(lv)
edgepoints=np.zeros(4,dtype=int)
ii=0

vold=0
eqpoifile = open('eqpointsrinzelcorI%d.txt' %Iapp, 'w')

for k in range(0,lv):
	vnew=w[np.where(abs(viso(w,v[k]))==min(abs(viso(w,v[k]))))[0][0]]
	if abs(viso(vnew,v[k]))<0.1:
		vvalues[k]=vnew
	else:
		vvalues[k]=vold
	if k>0:
		#if (winf(SV,v[k])-vold)*(winf(SV,v[k])-vnew)<0:
		#	eqpoifile.write('%.2f\n'%v[k])
		if (winf(SV,v[k])-vvalues[k])*(winf(SV,v[k-1])-vvalues[k-1])<=0:
			eqpoifile.write('%.2f\n'%v[k])
	if ii==0:
		if vold<0 and vnew>0 and v[k]<-10:
			edgepoints[ii]=k
			ii=ii+1
	if ii==1:
		if vold>0 and vnew<0 and v[k]>0:
			edgepoints[ii]=k+1
			ii=ii+1
	if ii==2:
		if vold<1 and vnew>1 and v[k]>10:
			edgepoints[ii]=k
			ii=ii+1
	if ii==3:
		if vold>0 and vnew<0 and v[k]>100:
			edgepoints[ii]=k	
	vold=vnew

eqpoifile.close() 
vrinzelfile = open('vncrinzel%d.txt' %Iapp, "w")
for kk in range(0,lv):
	vrinzelfile.write('%.4f %.4f\n'%(v[kk],vvalues[kk]))
vrinzelfile.close()
edgepointfile= open('edgerinzel%d.txt' %Iapp, "w")
for ll in range(0,4):
	edgepointfile.write('%d\n' %edgepoints[ll])
edgepointfile.close()

plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('recovery variable W')
plt.plot(v,winf(SV,v),label='W-nullcline')
plt.plot(v[edgepoints[0]:edgepoints[1]],vvalues[edgepoints[0]:edgepoints[1]],'orange',label='V-nullcline')
plt.plot(v[edgepoints[2]:edgepoints[3]],vvalues[edgepoints[2]:edgepoints[3]],'orange')
plt.ylim(0,1)
plt.legend()
plt.savefig('rinzelclinescorm15a62.pdf')

