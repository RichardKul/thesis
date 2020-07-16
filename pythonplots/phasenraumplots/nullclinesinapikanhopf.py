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

matplotlib.rcParams.update({'font.size': 22})
plt.axes().set_aspect(1)

style="Simple,tail_width=0.5,head_width=4,head_length=8"
kw = dict(arrowstyle=style, color="k")

a1 = patches.FancyArrowPatch((-36,0.87), (-70,0.5),connectionstyle="arc3,rad=.3",**kw )
a2 = patches.FancyArrowPatch((-70,0.45), (-55,0.15),connectionstyle="arc3,rad=.5",**kw)
a3 = patches.FancyArrowPatch((-52,0.15), (-30,0.7),connectionstyle="arc3,rad=.5", **kw)
a4 = patches.FancyArrowPatch((-30,0.75), (-35,0.85),connectionstyle="arc3,rad=.4", **kw)

I=46
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


t=np.arange(-75,-10,0.1)
plt.figure()

for a in [a1,a2,a3,a4]:
    plt.gca().add_patch(a)

plt.xlabel('membrane voltage V [mV]')
plt.ylabel('gating variable n')
plt.plot(t,vnc(t,I,vm,km,gL,EL,gNa,ENa,gK,EK),label='v-nullcline')
plt.plot(t,finf(t,vn,kn),label='n-nullcline')
plt.plot(-50.2138,finf(-50.2138,vn,kn),'bo')
plt.text(-50.2138-1,finf(-50.2138,vn,kn)+0.02,'$P$')
plt.xlim((-75,-10))
plt.ylim((-0.1,1.1))
plt.legend()
plt.tight_layout()
plt.savefig('inapikanhopfnc.pdf')

