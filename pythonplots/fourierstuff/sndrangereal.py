#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

date='realdrange6aem2'
nr=[9]
l=len(nr)
points=1000000
length=500000
dvalues=20
SNR=np.zeros((l,dvalues))
D=np.zeros((l,dvalues))
scale=np.zeros((l,dvalues))
ii=0

'''
for c in nr:
	for z in range(1,5):
		file=open('/home/richard/outhome/ft%s%d%d.txt' %(date,z,c),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])

for c in nr:
	for z in range(6,21):
		file=open('/home/richard/outhome/ft%s%d%d.txt' %(date,z,c),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		SNR[ii][z-2]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1
'''

for c in nr:
	for z in range(1,21):
		file=open('/home/richard/outhome/ft%s%d%d.txt' %(date,z,c),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date,z,c),"r")
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
		D[ii][z-1]=value[name.index('D')]
		dt=value[name.index('dt')]
		N=value[name.index('N')]-value[name.index('Neq')]
		repetitions=value[name.index('repetitions')]
		epsilon=value[name.index('epsilon')]
		omega=value[name.index('omega')]
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
plt.xlabel('noise intensity D')
plt.ylabel('SNR')
#xs=np.array([0.11,0.12,0.13,0.14,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.30])
#xs=np.arange(0.31,0.51,0.01)
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
for n in range(0,l):
	plt.plot(D[n,:],(SNR[n,:]-1)/scale[n,:],label='I=%f' %(nr[n]*0.02-0.1))
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
plt.legend()
plt.tight_layout()
#plt.plot(sax2,say2/T2,label='e6')
plt.savefig('snrauto%s.pdf' %date)
