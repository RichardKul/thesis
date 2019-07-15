#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/NetBeansProjects/realbursting/countI2.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax=np.array(x)
ay=np.array(y)
#az=40*np.array(z)-65
#aa=-np.array(a)

#area=(np.sum(ay)-min(ay)*points)*0.003/360

plt.figure()
plt.xlabel('time')
plt.ylabel('membrane voltage')
#plt.xlim(50,70)
plt.plot(ax,ay)
#plt.plot(ax,az)
#plt.plot(ax,aa)
plt.legend()
plt.savefig('realburst2.pdf')
