#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

import matplotlib.patches as patches

style="Simple,tail_width=0.5,head_width=4,head_length=8"
kw = dict(arrowstyle=style, color="k")

a1 = patches.FancyArrowPatch((5,1), (8,1),connectionstyle="arc3,rad=-.7",**kw )
a2 = patches.FancyArrowPatch((9.48,0.1), (3.2,0.1),connectionstyle="arc3,rad=0",**kw)
a3 = patches.FancyArrowPatch((3.2,0.2), (5,1),connectionstyle="arc3,rad=.3", **kw)
a4 = patches.FancyArrowPatch((8,1), (9.48,0.2),connectionstyle="arc3,rad=.3", **kw)


t=np.arange(2,10.5,0.001)
vnc=-5/2*np.sin(t)
xnc=0*t

plt.figure()
for a in [a1,a2,a3,a4]:
    plt.gca().add_patch(a)
plt.xlabel('position x')
plt.ylabel('velocity V')
plt.plot(t,xnc,label='x-Nullcline')
plt.plot(t,vnc,label='v-Nullcline')

plt.legend()
plt.savefig('nullclinemm2.pdf')

