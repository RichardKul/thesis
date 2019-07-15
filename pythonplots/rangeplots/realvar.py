#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

points=100000

file=open('/home/richard/NetBeansProjects/detmodelreal/voltageplusstate.txt',"r")
x,y,z,a=[],[],[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax=np.array(x)
ay=np.array(y)

v2=np.zeros(points)
vav=np.mean(ay)
for i in range(0,points):
	v2[i]=(ay[i]-vav)**2
v2av=np.mean(v2)	

plt.figure()
plt.xlabel('time')
plt.ylabel('membrane voltage')
plt.xlim(22000,23000)
plt.plot(ax,ay,label='%.6f' %v2av)
#plt.plot(ax,az)
#plt.plot(ax,aa)
plt.legend()
plt.savefig('realvar7jj.pdf')
