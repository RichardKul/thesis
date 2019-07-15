#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

date='realfrange9jj'
nr=[10]
l=len(nr)
omega0=0.0001
runs=160
points=1000000
length=500000
N=100000000
repetitions=20
dt=0.0001
T=N*repetitions*dt
S=np.zeros(length)

SNR=np.zeros((l,20))
ii=0

for c in nr:
	for z in range(1,21):
		file=open('/home/richard/outhome/ft%s%d%d.txt' %(date,c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		omega=omega0*z**2
		omegaind=round(omega*T)
		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
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
plt.xlabel('signal frequency')
plt.ylabel('SNR')
xs=omega0*np.square(np.arange(1,21))
#plt.yscale('log')
plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
for n in range(0,l):
	plt.plot(xs,SNR[n,:],label='I=0.08,D=%f' %(nr[n]*0.01))
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
plt.savefig('snrfrange2.pdf')
