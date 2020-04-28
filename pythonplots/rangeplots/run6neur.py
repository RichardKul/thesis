#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/mastergit/NetBeansProjects/realstatevar/realstate15.txt',"r")
x,y,z,a=[],[],[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
	z.append(float(row[3]))
	#a.append(float(row[4]))
ax=np.array(x)
ay=np.array(y)
az=40*np.array(z)-50
#aa=-np.array(a)/40

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
axs = plt.subplot(111)
plt.suptitle('$I$=$0.15\mu A/cm^2$')
#plt.suptitle('I=0')
plt.xlabel('time [s]')
plt.ylabel('membrane voltage [mV]')
plt.xlim(7.25,7.65)
axs.plot(ax/1000,ay,'black')
#plt.plot(ax/1000,az)
#plt.plot(ax,aa)
axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig('realstatevar15sh.pdf')
