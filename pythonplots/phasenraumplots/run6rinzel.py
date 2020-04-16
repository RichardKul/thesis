#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/mastergit/NetBeansProjects/realstatevar/realstaterinzelnb2.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)

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
Iapp=-15 #-12
SV=1.27

v=np.arange(-50,10,0.01)

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

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
axs = plt.subplot(111)
plt.suptitle('$I$=-15$\mu A/cm^2$')
#plt.suptitle('I=0')
#plt.xlabel('time [s]')
plt.ylabel('recovery variable W')
plt.xlabel('membrane voltage V [mV]')
plt.plot(colx[edge[0]:edge[1]],col[edge[0]:edge[1]],'blue',label='V-nullcline')
plt.plot(colx[edge[2]:edge[3]],col[edge[2]:edge[3]],'blue')
plt.plot(v,winf(SV,v),'orange',label='W-nullcline')
plt.xlim(-50,10)
axs.plot(ax,ay,'black')
plt.legend()
axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig('realstatedetrinzelpnb2wnblack.pdf')
