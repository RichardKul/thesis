#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt
axs=np.zeros(10000)
ays=np.zeros(10000)
azs=np.zeros(10000)
file=open('/home/richard/mainp2new/train29a202.txt',"r")
x,y,z=[],[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
	z.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)
az=np.array(z)
for z in range(0,10000):
	axs[z]=ax[z]
	ays[z]=ay[z]
	azs[z]=az[z]	

plt.figure()
plt.ylabel('gating variable n')
plt.xlabel('membrane voltage V')
plt.xlim(-26.1,-26)
plt.ylim(0.447,0.45)
plt.plot(ays,azs)
plt.savefig('mainp2foc.pdf')
