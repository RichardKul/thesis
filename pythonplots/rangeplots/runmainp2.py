#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt
axs=np.zeros(10000)
ays=np.zeros(10000)
azs=np.zeros(10000)
file=open('/home/richard/mainp2new/trainstatecomp5011.txt',"r")
x,y,z,a=[],[],[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
	z.append(float(row[3]))
	a.append(float(row[4]))
ax=np.array(x)
ay=np.array(y)
az=np.array(z)
aa=np.array(a)

#file2=open('/home/richard/mainp2new/train229a404.txt',"r")
#x2,y2,z2=[],[],[]
#for k in file2:
#	row=k.split()
#	x2.append(float(row[0]))
#	y2.append(float(row[1]))
#	z2.append(float(row[4]))
#ax2=np.array(x2)
#ay2=np.array(y2)
#az2=np.array(z2)
plt.figure()
plt.xlabel('time')
plt.ylabel('membrane voltage V')
#plt.xlim(850,950)
plt.plot(ax,ay,label='spike train')
plt.plot(ax,40*az-65,label='new bursting variable')
plt.plot(ax,40*aa-65,label='old bursting variable')
plt.legend()
plt.savefig('mainp2statecomp.pdf')
