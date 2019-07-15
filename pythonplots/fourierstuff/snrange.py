#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

omega=0.01
runs=500
points=1000000
length=500000
N=100000000
repetitions=2
dt=0.00001
T=N*repetitions*dt
T2=T*10
S=np.zeros(length)
omegaind=round(omega*T)
SNR=np.zeros((4,20))
ii=0

for c in [30,40,50]:
	for z in range(1,21):
		file=open('/home/richard/outhome/ft11je1%d%d.txt' %(c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1

repetitions=2
T=N*repetitions*dt

for c in [25]:
	for z in range(1,21):
		file=open('/home/richard/outhome/ft8je1%d%d.txt' %(c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind+1],ay[omegaind+2]])
	ii=ii+1
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
plt.xlabel('bias current')
plt.ylabel('SNR')
xs=np.arange(-0.75,4.25,0.25)
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
plt.plot(xs,SNR[0,:],label='D=3')
plt.plot(xs,SNR[1,:],label='D=4')
plt.plot(xs,SNR[2,:],label='D=5')
#plt.plot(xs,SNR[1,:],label='D=2.5')
plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
plt.savefig('snr345log.pdf')
