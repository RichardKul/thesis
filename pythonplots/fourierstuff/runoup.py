#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

epsilon=1
omega=10
tau=1
runs=1000
points=100000
length=50000
N=100000
dt=0.001
T=N*dt
T2=T*10
S=np.zeros(length)
amp=epsilon**2*T/(4*(1+(omega*tau)**2))

def spectrum(f,tau,sigma):
	return (2*tau*sigma**2)/(1+(2*np.pi*f*tau)**2)

file=open('/home/richard/NetBeansProjects/oup/xtraje5newwcav.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax=np.array(x)
ay=np.array(y)
ax2=np.zeros(length) 
ay2=np.zeros(length)
for l in range(0,length):
	ax2[l]=ax[l]
	ay2[l]=ay[l]

#files=open('/home/richard/NetBeansProjects/oup/xtraje6newwcav2.txt',"r")
#sx,sy=[],[]
#for ks in files:
#	rows=ks.split()
#	sx.append(float(rows[0]))
#	sy.append(float(rows[1]))
#sax=np.array(sx)
#say=np.array(sy)
#sax2=np.zeros(length) 
#say2=np.zeros(length)
#for l in range(0,length):
#	sax2[l]=sax[l]
#	say2[l]=say[l]
#plt.figure()
#plt.xlabel('time')
#plt.ylabel('position')
#plt.xlim(0,10)
#plt.plot(ax,ay)
#plt.plot(ax,aa)
#plt.savefig('oub.pdf')

#ys = fft(ay)
#for l in range(0,length):
#	S[l]=abs(ys[l])*abs(ys[l])/T

#omega=np.arange(0,length)*2*np.pi/T
plt.figure()
plt.xlabel('frequency')
plt.ylabel('power spectrum')
plt.yscale('log')
plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
plt.plot(ax2,ay2/T,label='sim')
plt.plot(ax2,spectrum(ax2,1,1),label='theory')
plt.plot(5/np.pi,amp,'o',label="delta")
#plt.plot(sax2,say2/T2,label='e6')
plt.legend()
plt.savefig('oupcopye5theordelta.pdf')
