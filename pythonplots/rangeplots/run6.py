#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/mastergit/NetBeansProjects/realstatevar/realstate155.txt',"r")
x,y,z,a=[],[],[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
	z.append(float(row[3]))
	a.append(float(row[4]))
ax=np.array(x)
ay=np.array(y)
az=40*np.array(z)-65
aa=-np.array(a)/40

plt.figure()
plt.suptitle('I=0.2$\mu A/cm^2$')
#plt.suptitle('I=0')
plt.xlabel('time [s]')
plt.ylabel('membrane voltage [mV]')
#plt.xlim(3.75,4.25)
plt.plot(ax/1000,ay,'black')
#plt.plot(ax,az)
#plt.plot(ax,aa)
plt.savefig('realstate155black.pdf')
