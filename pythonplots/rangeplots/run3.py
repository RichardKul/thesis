#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

points=100000

file=open('/home/richard/NetBeansProjects/inapik/inapppi0d0l.txt',"r")
x,y,z,a=[],[],[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax=np.array(x)
ay=np.array(y)
#az=40*np.array(z)-65
#aa=-np.array(a)

area=(np.sum(ay)-min(ay)*points)*0.003/360

plt.figure()
plt.xlabel('time')
plt.ylabel('membrane voltage')
plt.xlim(50,70)
plt.plot(ax,ay,label='%.6f' %area)
#plt.plot(ax,az)
#plt.plot(ax,aa)
plt.legend()
plt.savefig('inapburst8j.pdf')
