#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

def spectrum(f,tau,sigma):
	return (2*tau*sigma**2)/(1+(2*np.pi*f*tau)**2)

omega=0.000005
runs=50
points=1000000
length=500000

date='realrinzelpoi26n1'

D=500
istart=1
ivalues=4
nr=9
S=np.zeros(length)
#omegaind=round(omega*T)
for z in range(istart,istart+ivalues):
	file=open('/home/richard/outhome/spike%s%d%d.txt' %(date,D,z),"r")
	#file=open('/home/richard/outhome/spikerealrinzel25o%d%d.txt' %(z,nr),"r")
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
	param=open('/home/richard/outhome/param%s%d%d.txt' %(date,D,z),"r")
	ll=0
	name,value=[],[]
	for k in param:
		row=k.split()
		lp=len(row)
		if ll<1:
			for jj in range(0,lp):
				name.append(row[jj])
		else:
			for kk in range(0,lp):
				value.append(float(row[kk]))
		ll=ll+1
	dt=value[name.index('dt')]
	N=value[name.index('N')]-value[name.index('Neq')]
	repetitions=value[name.index('repetitions')]
	T=N*repetitions*dt

#background=np.mean([ay2[omegaind-1],ay2[omegaind-2],ay2[omegaind+1],ay2[omegaind+2],ay2[omegaind-3]])
#SNR=ay2[omegaind]/background
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
	ax3=np.arange(0,length)
#omega=np.arange(0,length)*2*np.pi/T
	plt.figure()
	#plt.suptitle('I=%.2f, D=%.2f' %(-0.1+z*0.02,D*0.01))
	plt.suptitle('I=%.2f, D=%.2f' %(-12+z*0.6,D/10))
	plt.xlabel('Frequency $[10^3s^{-1}]$')
	plt.ylabel('Spectral power')	
	plt.yscale('log')
	plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
	plt.plot(ax2,ay2/T,label='Simulation')
#plt.plot(ax2,spectrum(ax2,1/(4.748*10**(-3)),13),label='theory')
#plt.plot(omega,background/T,'kx')
	#plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
	#plt.savefig('inapikrealfrange9aspD=%.2fI=%.2f.pdf' %(D*0.01,-0.1+z*0.02))
	plt.savefig('inapik%sD=%.2fI=%.2f.pdf' %(date,D/10,-12+z*0.6))	
