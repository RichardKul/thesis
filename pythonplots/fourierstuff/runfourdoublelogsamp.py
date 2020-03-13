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

timefac=1000

date='realanhopf25fsamp'

mode='spike'
D=35
I=10
istart=0
ivalues=10
nr=9
#zv=[20,30,40,50,60,70,80,90,100,200]
zv=np.arange(1,11)

for z in range(istart,istart+ivalues):
	freq=200*zv[z]
	#freq=200*z
	file=open('/home/richard/outhome/%s%s%d%d%d.txt' %(mode,date,D,I,zv[z]),"r")
	#file=open('/home/richard/outhome/%s%s%d%d%d.txt' %(mode,date,D,I,z),"r")
	#file=open('/home/richard/outhome/spikerealrinzel25o%d%d.txt' %(z,nr),"r")
	x,y=[],[]
	for k in file:
		row=k.split()
		x.append(float(row[0]))
		y.append(float(row[1]))
	ax=np.array(x)
	ay=np.array(y)
	l=len(ax)
	if l==0:
		continue
	#file2=open('/home/richard/outhome/%sshort%s%d%d%d.txt' %(mode,date,D,I,z),"r")
	file2=open('/home/richard/outhome/%sshort%s%d%d%d.txt' %(mode,date,D,I,zv[z]),"r")
	#file=open('/home/richard/outhome/spikerealrinzel25o%d%d.txt' %(z,nr),"r")
	x,y=[],[]
	for k in file2:
		row=k.split()
		x.append(float(row[0]))
		y.append(float(row[1]))
	ax2=np.array(x)
	ay2=np.array(y)
	l2=len(ax2)
	if l2==0:
		continue
	#param=open('/home/richard/outhome/param%s%d%d%d.txt' %(date,D,I,z),"r")
	param=open('/home/richard/outhome/param%s%d%d%d.txt' %(date,D,I,zv[z]),"r")
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
	#rate=open('/home/richard/outhome/g%s%d%d%d.txt' %(date,D,I,z),"r")
	rate=open('/home/richard/outhome/g%s%d%d%d.txt' %(date,D,I,zv[z]),"r")
	for k in rate:
		row=k.split()
		r=float(row[1])
	rate0=open('/home/richard/outhome/grealanhopf22j3020.txt',"r")
	for k in rate0:
		row=k.split()
		r0=float(row[1])
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
	ind11=np.argmin(abs(ax-r0))
	ind21=np.argmin(abs(ax2-r0))
	ind1=np.argmax(ay[ind11-5:ind11+5])+ind11-5	
	ind2=np.argmax(ay2[ind21-5:ind21+5])+ind21-5
#ys = fft(ay)
#for l in range(0,length):
#	S[l]=abs(ys[l])*abs(ys[l])/T
#omega=np.arange(0,length)*2*np.pi/T
	plt.figure()
	#plt.suptitle('I=%.2f, D=%.2f' %(-0.1+z*0.02,D*0.01))
	plt.suptitle('I=%.2f, D=%.2f' %(45.5,D/100))
	plt.xlabel('Frequency $[s^{-1}]$')
	plt.ylabel('Spectral power')	
	plt.yscale('log')
	plt.xscale('log')
#	plt.xlim(ax[0]*timefac,10000000/(2*T)*timefac)
#plt.xlim(4*10**(-4),100)
	plt.plot(ax[0:l-11]*timefac,ay[0:l-11]*ax[0]*timefac,'blue')#,label='Simulation')
	plt.plot(ax[l-11:l]*timefac,ay[l-11:l]*ax[0]*timefac,'blue')
	plt.plot(ax2*timefac,ay2*ax2[0]*timefac,'orange')#,label='Simulation')
	plt.plot(ax[ind1]*timefac,ay[ind1]*timefac*ax[0],'x',color='blue')
	plt.plot(ax2[ind2]*timefac,ay2[ind2]*timefac*ax2[0],'x',color='orange')
	plt.plot(ax[0:l-11]*timefac,r*np.ones(l-11)*timefac,'g',label='running firing rate')
	plt.plot(ax2[0:l2]*timefac,r*np.ones(l2)*timefac,'g')
#plt.plot(ax2,spectrum(ax2,1/(4.748*10**(-3)),13),label='theory')
#plt.plot(omega,background/T,'kx')
	#plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
	#plt.savefig('inapikrealfrange9aspD=%.2fI=%.2f.pdf' %(D*0.01,-0.1+z*0.02))
	plt.savefig('inapikanhopf%s2%sfourierD=%.2fI=%.2f.pdf' %(mode,date,D/100,freq))	
