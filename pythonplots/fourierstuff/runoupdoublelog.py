#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft



def spectrum(f,tau,sigma):
	return (2*tau*sigma**2)/(1+(2*np.pi*f*tau)**2)

date='xtraje18ftfteneq2'


nr=9
epsilon=0.01
omega=0.01
tau=3



file=open('/home/richard/mastergit/NetBeansProjects/oup/%s.txt' %(date),"r")
#file=open('/home/richard/outhome/spikerealrinzel25o%d%d.txt' %(z,nr),"r")
x,y=[],[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax=np.array(x)
ay=np.array(y)
file2=open('/home/richard/mastergit/NetBeansProjects/oup/%sshort.txt' %(date),"r")
#file=open('/home/richard/outhome/spikerealrinzel25o%d%d.txt' %(z,nr),"r")
x,y=[],[]
for k in file2:
	row=k.split()
	x.append(float(row[0]))
	y.append(float(row[1]))
ax2=np.array(x)
ay2=np.array(y)

amp=epsilon**2/ax[0]/(4*(1+((omega/2)*tau)**2))
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
#omega=np.arange(0,length)*2*np.pi/T
plt.figure()
plt.xlabel('Frequency $[s^{-1}]$')
plt.ylabel('Spectral power')	
plt.yscale('log')
plt.xscale('log')
plt.xlim(10**(-4),950)
plt.ylim(10**(-10),10**3)
#	plt.xlim(ax[0]*timefac,10000000/(2*T)*timefac)
#plt.xlim(4*10**(-4),100)
plt.plot(ax,ay*ax[0],label='Simulationlong')
plt.plot(ax2,ay2*ax2[0],label='Simulationshort')
plt.plot(ax,spectrum(ax,3,5),label='theorylong')
plt.plot(ax2,spectrum(ax2,3,5),label='theoryshort')
plt.plot(omega/2,amp+spectrum(omega/(4*np.pi),3,5),'x',label="delta")
plt.legend()
#plt.plot(ax2,spectrum(ax2,1/(4.748*10**(-3)),13),label='theory')
#plt.plot(omega,background/T,'kx')
	#plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
	#plt.savefig('inapikrealfrange9aspD=%.2fI=%.2f.pdf' %(D*0.01,-0.1+z*0.02))
plt.savefig('oup%sfourier.pdf' %(date))	
