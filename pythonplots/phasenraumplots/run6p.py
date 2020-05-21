#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/mastergit/NetBeansProjects/realstatevar/realstate125.txt',"r")
x,y,z,a=[],[],[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
	z.append(float(row[3]))
	a.append(float(row[4]))
ax=np.array(x)
ay=np.array(y)
az=40*np.array(z)-65
aa=-np.array(a)/40

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
t=np.arange(-80,0,0.01)

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
axs = plt.subplot(111)
plt.plot(t,vnc(t,I,vm,km,gL,EL,gNa,ENa,gK,EK),label='v-nullcline')
plt.plot(t,finf(t,vn,kn),label='n-nullcline')
plt.suptitle('$I$=-0.1$\mu A/cm^2$')
plt.ylabel('gating variable n')
plt.xlabel('membrane voltage [mV]')
#plt.xlim(3.75,4.25)
axs.plot(ax[25001:28433],ay[25001:28433],'black')
#plt.plot(ax,az)
#plt.plot(ax,aa)
axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
plt.legend()
plt.tight_layout()
plt.savefig('realstatep125vshbig.pdf')
