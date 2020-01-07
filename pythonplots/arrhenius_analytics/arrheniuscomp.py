#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

date0='raterealrinzelpoi26n1realrinzel15ninv00'
file0=open('%s.txt'%date0,"r")
x0=[]
for k in file0:
	row=k.split()
	x0.append(float(row[0]))
ax0=np.array(x0)

datem1='raterealrinzel25onewrealfast19jjem2st6'
filem1=open('%s.txt'%datem1,"r")
xm1=[]
for k in filem1:
	row=k.split()
	xm1.append(float(row[0]))
axm1=np.array(xm1)

datep1='raterealrinzel25onewrealfast19jjem2st9'
filep1=open('%s.txt'%datep1,"r")
xp1=[]
for k in filep1:
	row=k.split()
	xp1.append(float(row[0]))
axp1=np.array(xp1)

y=np.array([ax0[1],ax0[3],ax0[5]])
y2=np.array([ax0[7],ax0[9],ax0[11]])
x=np.array([1/20,1/30,1/50])
plt.figure()
plt.ylabel('transition rate')
plt.xlabel('inverse noise intensity')
plt.yscale('log')
plt.plot(x,1/y,label='btoeq I=-11.4 with signal')
plt.plot(x,1/axm1[0:3],label='btoeq I=-12.4 without signal')
plt.plot(x,1/axp1[0:3],label='btoeq I=-10 without signal')
plt.plot(x,1/y2,label='eqtob I=-11.4 with signal')
plt.plot(x,1/axm1[3:6],label='eqtob I=-12.4 without signal')
plt.plot(x,1/axp1[3:6],label='eqtob I=-10 without signal')
plt.legend()
plt.savefig('arrhcomp.pdf')

