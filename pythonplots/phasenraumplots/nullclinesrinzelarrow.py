#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches

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

plt.axes().set_aspect(1)

style="Simple,tail_width=0.5,head_width=4,head_length=8"
kw = dict(arrowstyle=style, color="k")

a1 = patches.FancyArrowPatch((55,1.07), (13,1),connectionstyle="arc3,rad=.3",**kw )
a2 = patches.FancyArrowPatch((3,0.9), (8,0.7),connectionstyle="arc3,rad=.2",**kw)
a3 = patches.FancyArrowPatch((15,0.7), (75,0.85),connectionstyle="arc3,rad=.3", **kw)
a4 = patches.FancyArrowPatch((75,0.9), (60,1.05),connectionstyle="arc3,rad=.4", **kw)



v=np.arange(-50,150,0.01)

plt.figure()

for a in [a1,a2,a3,a4]:
    plt.gca().add_patch(a)

plt.xlabel('membrane voltage V [mV]')
plt.ylabel('recovery variable W')
plt.plot(v,winf(SV,v),'orange',label='W-nullcline')
plt.plot(colx[edge[0]:edge[1]],col[edge[0]:edge[1]],'blue',label='V-nullcline')
plt.plot(colx[edge[2]:edge[3]],col[edge[2]:edge[3]],'blue')
plt.plot(-23.2,winf(SV,-23.2),'bo')
plt.plot(1.3,winf(SV,1.3),'o',color='purple')
plt.plot(20.9,winf(SV,20.9),'ro')
plt.text(-23.2-7,winf(SV,-23.2)+0.02,'$P_1$')
plt.text(1.3-6,winf(SV,1.3)+0.02,'$P_2$')
plt.text(20.9-5,winf(SV,20.9)+0.02,'$P_3$')
plt.ylim(0,1.1)
plt.legend()
plt.savefig('rinzelclinesarrowwp.pdf')

