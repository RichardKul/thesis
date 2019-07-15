#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

file=open('/home/richard/NetBeansProjects/cpphellodet/mechdet2.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[1]))
	y.append(float(row[2]))
ax=np.array(x)
ay=np.array(y)

plt.figure()
plt.ylabel('velocity v')
plt.xlabel('position x')
plt.plot(ay,ax)
plt.savefig('mechdet2.pdf')
