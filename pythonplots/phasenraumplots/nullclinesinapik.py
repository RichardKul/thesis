#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def finf(x,v12,k):
	return 1/(1+np.exp((v12-x)/k))
def vnc(x,I,v12,k,gL,EL,gNa,ENa,gK,EK):
	return (I-gL*(x-EL)-gNa*finf(x,v12,k)*(x-ENa))/(gK*(x-EK)) 

plt.axes().set_aspect(1)

style="Simple,tail_width=0.5,head_width=4,head_length=8"
kw = dict(arrowstyle=style, color="k")

a1 = patches.FancyArrowPatch((-17,0.85), (-40,0.3),connectionstyle="arc3,rad=.2",**kw )
a2 = patches.FancyArrowPatch((-40,0.23), (-36,0.1),connectionstyle="arc3,rad=.5",**kw)
a3 = patches.FancyArrowPatch((-34,0.1), (-11,0.72),connectionstyle="arc3,rad=.5", **kw)
a4 = patches.FancyArrowPatch((-11,0.75), (-15,0.85),connectionstyle="arc3,rad=.4", **kw)

I=0
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


t=np.arange(-75,-10,0.1)
plt.figure()

for a in [a1,a2,a3,a4]:
    plt.gca().add_patch(a)

plt.xlabel('membrane voltage V [mV]')
plt.ylabel('gating variable n')
plt.plot(t,vnc(t,I,vm,km,gL,EL,gNa,ENa,gK,EK),label='v-nullcline')
plt.plot(t,finf(t,vn,kn),label='n-nullcline')
plt.plot(-69.108,finf(-69.108,vn,kn),'bo')
plt.plot(-55.8294,finf(-55.8294,vn,kn),'o',color='purple')
plt.plot(-21.7225,finf(-21.7225,vn,kn),'ro')
plt.text(-69.108,finf(-69.108,vn,kn)+0.02,'$P_1$')
plt.text(-55.8294-2,finf(-55.8294,vn,kn)+0.02,'$P_2$')
plt.text(-21.7225-3,finf(-21.7225,vn,kn),'$P_3$')
plt.xlim((-75,-10))
plt.ylim((-0.1,0.9))
plt.legend()
plt.savefig('inapikrealncwnp.pdf')

