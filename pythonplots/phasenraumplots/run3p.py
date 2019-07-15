#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/NetBeansProjects/inapik/inapi0d0.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)

plt.figure()
plt.xlabel('membrane voltage V')
plt.ylabel('gating variable n')
#plt.gcf().subplots_adjust(left=0.15)
plt.plot(ax,ay)
plt.savefig('inapburst.pdf')
