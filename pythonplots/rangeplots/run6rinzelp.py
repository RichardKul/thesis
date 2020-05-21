#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/mastergit/NetBeansProjects/realstatevar/realstaterinzel16100.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax=np.array(x)
ay=np.array(y)

matplotlib.rcParams.update({'font.size': 22})

plt.figure()
axs = plt.subplot(111)
plt.suptitle('$I$=-16$\mu A/cm^2$')
#plt.suptitle('I=0')
#plt.xlabel('time [s]')
plt.ylabel('recovery variable W')
plt.xlabel('membrane voltage V [mV]')
#plt.xlim(0,0.1)
axs.plot(ax/1000,ay)
axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)
plt.tight_layout()
plt.savefig('realstatedetrinzel16200ps.pdf')
