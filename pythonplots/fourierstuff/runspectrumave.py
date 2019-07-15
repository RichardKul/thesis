#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import cmath

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

ave=10
dt=0.00001
omega=0.001
runs=10
repetitions=50
fpoints=50000
datapoints=1000000
N=100000000
T=N*repetitions*dt
S=np.zeros(fpoints)
last_entry=0
def exp(f,t):
	return cmath.exp(2*np.pi*1j*f*t)

file=open('/home/richard/outhome/spike8je1long3010.txt',"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
ax=np.array(x)
for z in range(0,datapoints):
	current=ax[z]
	if current>last_entry:
		last_entry=current
	else:
		points_act=z
		break

length=int(points_act/ave)
x_act=np.zeros(length)
for hh in range(0,ave):
	for z2 in range(0,length):
		x_act[z2]=ax[z2+hh*length]
	for ii in range(0,fpoints):
		f=3*ii/fpoints
		arr=np.zeros(length,dtype=np.complex_)
		for jj in range(0,length):
			arr[jj]=exp(f,ax[jj])
		S[ii]=(hh*S[ii]+(np.abs(np.sum(arr))**2)/T)/(hh+1)	

f=np.arange(0,3,3/fpoints)
plt.figure()
plt.xlabel('frequency')
plt.ylabel('power spectrum')
plt.yscale('log')
plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
plt.plot(f,S)
#plt.plot(sax2,say2/T2,label='e6')
plt.savefig('deltaspectrum2.pdf')
