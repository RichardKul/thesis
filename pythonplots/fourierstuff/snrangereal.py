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

date='realfast9aem2sh'
date1='realfast19jjem2st'
D=[35]
D1=[]
Dtot=D+D1
l=len(D)+len(D1)
omega=0.000005
runs=50
points=1000000
length=500000
N=200000000
repetitions=200
dt=0.0005
epsilon=0.01
T=N*repetitions*dt
scale=T*epsilon
T1=T
T2=T*10
S=np.zeros(length)
omegaind=round(omega*T)
omegaind1=round(omega*T1)
SNR=np.zeros((l,20))
ii=0

for c in D:
	for z in range(1,21):
		file=open('/home/richard/outhome/ft%s%d%d.txt' %(date,c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1

for c1 in D1:
	for z in range(1,21):
		file=open('/home/richard/outhome/ft%s%d%d.txt' %(date1,c1,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		SNR[ii][z-1]=ay[omegaind1]/np.mean([ay[omegaind1-1],ay[omegaind1-2],ay[omegaind1+1],ay[omegaind1+2],ay[omegaind1-3]])
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

t=np.arange(-0.1,0.3,0.01)
xs=np.arange(-0.08,0.32,0.02)
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
for n in range(0,l):
	plt.plot(xs,(SNR[n,:]-1)/scale,label='D=%s' %(Dtot[n]*0.01))
for n in range(0,l):	
	plt.plot(t,400*snr(r0p,r0m,ap,am,bp,bm,t,Dtot[n]*0.01))
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [0,2,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
plt.legend()
plt.savefig('snrreal9awf.pdf')
