#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.xlim(-5,5)
plt.ylim(-9,7)
plt.axes().set_aspect(1)

style="Simple,tail_width=0.5,head_width=4,head_length=8"
kw = dict(arrowstyle=style, color="k")

a1 = patches.FancyArrowPatch((-4,-6), (0,6),**kw )
a2 = patches.FancyArrowPatch((0,6), (4,-6),**kw)
a3 = patches.FancyArrowPatch((-4,-6), (4,-6),connectionstyle="arc3,rad=.5", **kw)

for a in [a1,a2,a3]:
    plt.gca().add_patch(a)
plt.savefig('test2.pdf')
