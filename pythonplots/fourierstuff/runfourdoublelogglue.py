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

date='realfast15mtf'

mode='spike'
D=45
istart=6
ivalues=1
nr=9

for z in range(istart,istart+ivalues):
	file=open('/home/richard/outhome/%s%s%d%d.txt' %(mode,date,D,z),"r")
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
	file2=open('/home/richard/outhome/%sshort%s%d%d.txt' %(mode,date,D,z),"r")
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
	rate=open('/home/richard/outhome/g%s%d%d.txt' %(date,D,z),"r")
	for k in rate:
		row=k.split()
		r=float(row[1])
	rate0=open('/home/richard/outhome/grealfast3jje02020.txt',"r")
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
	axmax=np.argmin(abs(ax-ax2[10]))+1
	ind11=np.argmin(abs(ax-r0))
	ind21=np.argmin(abs(ax2-r0))
	ind1=np.argmax(ay[ind11-5:ind11+5])+ind11-5	
	ind2=np.argmax(ay2[ind21-5:ind21+5])+ind21-5
#ys = fft(ay)
#for l in range(0,length):
#	S[l]=abs(ys[l])*abs(ys[l])/T
#omega=np.arange(0,length)*2*np.pi/T
	plt.figure()
	plot1 = plt.subplot(111)
	plt.suptitle('I=%.2f, D=%.2f' %(-0.1+z*0.02,D*0.01))
	#plt.suptitle('I=%.2f, D=%.2f' %(43+z*0.25,D/100))
	plt.xlabel('Frequency $[s^{-1}]$')
	plt.ylabel('Spectral power $[s^{-1}]$')	
	plt.yscale('log')
	plt.xscale('log')
#	plt.xlim(ax[0]*timefac,10000000/(2*T)*timefac)
#plt.xlim(4*10**(-4),100)
	plt.plot(ax[0:axmax]*timefac,ay[0:axmax]*ax[0]*timefac,'blue')#,label='Simulation')
	plt.plot(ax[l-11:l]*timefac,ay[l-11:l]*ax[0]*timefac,'blue')
	plt.plot(ax2[10:-1]*timefac,ay2[10:-1]*ax2[0]*timefac,'blue')#,label='Simulation')
	#plt.plot(ax[ind1]*timefac,ay[ind1]*timefac*ax[0],'x',color='blue')
	#plt.plot(ax2[ind2]*timefac,ay2[ind2]*timefac*ax2[0],'x',color='orange')
	plt.plot(ax[0:axmax]*timefac,r*np.ones(axmax)*timefac,'g',label='overall firing rate')
	plt.plot(ax2[0:l2]*timefac,r*np.ones(l2)*timefac,'g')
	plt.plot(ax[0:axmax]*timefac,5200*np.ones(axmax),'r--',label='background')
	plt.plot(ax2[0:l2]*timefac,5200*np.ones(l2),'r--')
	#plt.plot(ax[0:axmax]*timefac,3*np.ones(axmax),'r--',label='background')
	#plt.plot(ax2[0:l2]*timefac,3*np.ones(l2),'r--')
	#plt.arrow(5*10**(-3),300,0,-299,length_includes_head=True)
	plt.arrow(0.222, 0.75, 0, -0.75, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	#plt.arrow(0.222, 0.23, 0, -0.23, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	#plt.text(10**(-2),1,'signal')
	#plt.text(10**(-2),0.7,'frequency')
	plt.text(10**(-2),0.5,'signal')
	plt.text(10**(-2),0.2,'frequency')
	#plt.arrow(70,10,0,-9,length_includes_head=True)
	#plt.arrow(0.5862, 0.6, 0, -0.6, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	plt.arrow(0.5862, 0.3, 0, -0.3, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	#plt.text(10**2,1,'firing rate')
	#plt.text(10**(2),0.7,'in spiking state')
	plt.text(10**2,0.5,'firing rate')
	plt.text(10**(2),0.2,'in spiking state')
#plt.plot(ax2,spectrum(ax2,1/(4.748*10**(-3)),13),label='theory')
#plt.plot(omega,background/T,'kx')
	plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
	#plt.savefig('specglue%s%sD=%.2fI=%.2f.pdf' %(mode,date,D*0.01,-0.1+z*0.02))
	plt.savefig('specpaper3.pdf')
	#plt.savefig('inapikanhopf%s2%sfourierD=%.2fI=%.2f.pdf' %(mode,date,D/100,43+z*0.25))	
