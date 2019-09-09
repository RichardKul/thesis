#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

date='realfrange8aem2'
nr=[30]
l=len(nr)

points=1000000
length=500000

fvalues=20
SNR=np.zeros((l,fvalues))
omegavec=np.zeros((l,fvalues))
scale=np.zeros((l,fvalues))
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
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date,c,z),"r")
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
		epsilon=value[name.index('epsilon')]
		omega=value[name.index('omega')]
		omegavec[ii][z-1]=omega
		T=N*repetitions*dt
		scale[ii][z-1]=epsilon**2*T
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
#xs=omega0*np.square(np.arange(1,21))
plt.yscale('log')
plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.ylim(1,100)
for n in range(0,l):
	plt.plot(omegavec[n,:],(SNR[n,:]-1)/scale[n,:],label='I=0.08,D=%f' %(nr[n]*0.01+0.1))
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
plt.savefig('snrautofrange8a.pdf')
