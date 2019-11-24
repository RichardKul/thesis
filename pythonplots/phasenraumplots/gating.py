#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def finf(t,v12,k12):
	return 1/(1+np.exp((v12-t)/k12))


t=np.arange(-70,30,0.01)

plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable')
plt.plot(t,finf(t,-18,14),label='$m_\infty(V)$')
plt.plot(t,finf(t,-25,5),label='$n_\infty(V)$')
plt.legend()
plt.savefig('gating.pdf')

