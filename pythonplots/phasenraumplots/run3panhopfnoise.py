#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def finf(x,v12,k):
	return 1/(1+np.exp((v12-x)/k))
def vnc(x,I,v12,k,gL,EL,gNa,ENa,gK,EK):
	return (I-gL*(x-EL)-gNa*finf(x,v12,k)*(x-ENa))/(gK*(x-EK)) 

I=45
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

file=open('/home/richard/mastergit/NetBeansProjects/realstatevar/realstateanhopf5.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)
#ayi=ay[::-1]

t=np.arange(-70,0,0.1)

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
axs = plt.subplot(111)
plt.suptitle('$I$=45$\mu A/cm^2$')
#plt.xlabel('time [s]')
plt.xlabel('membrane voltage V [mV]')
plt.ylabel('gating variable n')
plt.plot(t,vnc(t,I,vm,km,gL,EL,gNa,ENa,gK,EK),label='v-nullcline')
plt.plot(t,finf(t,vn,kn),label='n-nullcline')
#plt.gcf().subplots_adjust(left=0.15)
#plt.plot(abs(ax[0:1000])/1000,ayi[0:1000],'black')
#plt.xlim(39.97,40.17)
axs.plot(ax[4*33770:4*33970],ay[4*33770:4*33970],'black')
#plt.xlim(0,0.2)
#plt.xlim(0,0.05)
#plt.ylim(0.1,0.4)
plt.legend()
axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig('inaprealanhopfreal5psh.pdf')
