#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

t=np.arange(-68.44,-68.4,0.001)
vnc=((-1-(t+70))*(1+np.exp((-40-t)/15))**3)/(10*(t-60))
nnc=1/(1+np.exp((42+t)/12))
plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
plt.plot(t,vnc,label='v-nullcline')
plt.plot(t,nnc,label='n-nullcline')
plt.legend()
plt.savefig('inatnc1.pdf')

