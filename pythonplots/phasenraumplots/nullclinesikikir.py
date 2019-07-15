#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

t=np.arange(-70,-30,0.001)
vnc=65/(2*(t+80))-20/((1+np.exp((-80-t)/(-12)))*2)
nnc=1/(1+np.exp((-40-t)/5))
plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
plt.plot(t,vnc,label='v-nullcline')
plt.plot(t,nnc,label='n-nullcline')
plt.legend()
plt.savefig('ikikirnc.pdf')

