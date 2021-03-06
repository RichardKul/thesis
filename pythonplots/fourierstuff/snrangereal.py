#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def snr(r0p,r0m,ap,am,bp,bm,t,D):
	return (ap-am)**2/(D**2*(r(r0p,ap,bp,t,D)**(-1)+r(r0m,am,bm,t,D)**(-1)))

bp = 1.74
ap = 5.64
r0p = 0.0075
bm = 3.15
am = -10.76
r0m = 0.012

date='realrinzel20ninv0'
date2='realrinzel20ninv1'
date1='realrinzel15ninv0'
date3='realrinzel15ninv1'
D=[200,300]
D1=[500]
D2=[200,300]
D3=[500]
Dtot=D+D1+D2+D3
l=len(D)+len(D1)+len(D2)+len(D3)

points=1000000
length=500000
istart=2
ivalues=8
SNR=np.zeros((l,ivalues))
scale=np.zeros((l,ivalues))
ii=0

for c in D:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date,c,z),"r")
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
		T=N*repetitions*dt
		scale[ii][z-istart]=epsilon**2*T
		omegaind=round(omega*T)
		if c==200 and z==9:
			omegaind=2105		
		SNR[ii][z-istart]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1

for c1 in D1:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date1,c1,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date1,c1,z),"r")
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
		T=N*repetitions*dt
		scale[ii][z-istart]=epsilon**2*T
		omegaind=round(omega*T)		
		SNR[ii][z-istart]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1

for c2 in D2:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date2,c2,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		if len(x)==0:
			SNR[ii][z-istart]=SNR[ii][z-istart-1]
			scale[ii][z-istart]=scale[ii][z-istart-1]
		else:
			param=open('/home/richard/outhome/param%s%d%d.txt' %(date2,c2,z),"r")
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
			T=N*repetitions*dt
			scale[ii][z-istart]=epsilon**2*T
			omegaind=round(omega*T)		
			SNR[ii][z-istart]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1

for c3 in D3:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date3,c3,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date3,c3,z),"r")
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
		T=N*repetitions*dt
		scale[ii][z-istart]=epsilon**2*T
		omegaind=round(omega*T)		
		SNR[ii][z-istart]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
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
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('SNR')

t=np.arange(-0.1,0.3,0.01)
xs=np.arange(-15.4,-9,0.8)
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
plt.ylim(10**(-7),10**(-3))
colorv=['y','g','b','r','c']
for n in range(0,3):
	plt.plot(xs,abs((SNR[n,:]-1)/scale[n,:]),label='D=%.2f, rest' %(Dtot[n]*0.1))
for n in range(3,l):
	plt.plot(xs,abs((SNR[n,:]-1)/scale[n,:]),label='D=%.2f, run' %(Dtot[n]*0.1))
#for n in range(0,l):	
#	plt.plot(t,snr(r0p,r0m,ap,am,bp,bm,t,Dtot[n]*0.01),colorv[n]+'o')
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [0,2,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
plt.legend()
plt.savefig('snronlyrinzel%s.pdf' %date)
