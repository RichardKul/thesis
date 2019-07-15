#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


#plt.axes().set_aspect(1)

#style="Simple,tail_width=0.5,head_width=4,head_length=8"
#kw = dict(arrowstyle=style, color="k")

#a1 = patches.FancyArrowPatch((-24,0.6), (-38,0.33),connectionstyle="arc3,rad=.1",**kw )
#a2 = patches.FancyArrowPatch((-40,0.2), (-36,0.1),connectionstyle="arc3,rad=.5",**kw)
#a3 = patches.FancyArrowPatch((-34,0.05), (-20,0.45),connectionstyle="arc3,rad=.5", **kw)
#a4 = patches.FancyArrowPatch((-20,0.5), (-23,0.6),connectionstyle="arc3,rad=.5", **kw)

t=np.arange(-70,-17,0.1)
vnc=(0-8*(t+80)-20*(t-60)/(1+np.exp((-20-t)/15)))/(9*(t+90))
nnc=1/(1+np.exp((-25-t)/5))

plt.figure()

#for a in [a1,a2,a3,a4]:
#    plt.gca().add_patch(a)

plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
plt.plot(t,vnc,label='v-nullcline')
plt.plot(t,nnc,label='n-nullcline')
plt.xlim((-70,-17))
plt.ylim((-0.1,0.8))
plt.legend()
plt.savefig('inapikncpure.pdf')

