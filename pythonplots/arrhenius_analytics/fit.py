#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


def fit(x,x0,sigma):
	return np.exp(-(x-x0)**2/(2*sigma**2))/((2*sigma**2)**(0.5))

x=np.arange(-0.2,0.3,0.01)
plt.plot(x,fit(x,0.08,0.09))
plt.yscale('log')
plt.savefig('fit.pdf')
