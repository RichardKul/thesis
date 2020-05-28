#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

def spectrum(f,tau,sigma):
	return (2*tau*sigma**2)/(1+(2*np.pi*f*tau)**2)

matplotlib.rcParams.update({'font.size': 16})

omega=0.000005
runs=50

timefac=1000

date='realrinzelrange22mtf1'

mode='spike'
D=500
istart=8
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
	deff=open('/home/richard/outhome/d%s%d%d.txt' %(date,D,z),"r")
	for k in deff:
		row=k.split()
		d=float(row[1])
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
	plot1= plt.subplot(111)
	plt.suptitle('I=%.2f, D=%.0f' %(-17+z*0.8,D*0.1))
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
	#plt.plot(ax[0:axmax]*timefac,18*np.ones(axmax),'r--',label='background')
	#plt.plot(ax2[0:l2]*timefac,18*np.ones(l2),'r--')
	#plt.plot(ax[0:axmax]*timefac,2*d*np.ones(axmax)*timefac,'y',label='$2D_{eff}$')
	#plt.plot(ax2[0:l2]*timefac,2*d*np.ones(l2)*timefac,'y')
	#xdisplay, ydisplay = plot1.transAxes.transform_point((1, 2*d*timefac))
	#plt.arrow(0.1, 0.1, -0.1, 0,transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	#plt.plot(ax[0:axmax]*timefac,30000*np.ones(axmax),'r--',label='background')
	#plt.plot(ax2[0:l2]*timefac,30000*np.ones(l2),'r--')
	plt.plot(ax[0:axmax]*timefac,170*np.ones(axmax),'r--',label='background')
	plt.plot(ax2[0:l2]*timefac,170*np.ones(l2),'r--')
	#plt.arrow(5*10**(-3),300,0,-299,length_includes_head=True)
	#plt.arrow(0.277, 0.83, 0, -0.83, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	#plt.arrow(0.3457, 0.55, 0, -0.55, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	plt.arrow(0.3457, 0.26, 0, -0.26, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	plt.text(3*10**(-4),40,'signal')
	plt.text(3*10**(-5),20,'frequency')
	#plt.text(4*10**(-5),0.1,'signal')
	#plt.text(4*10**(-5),0.04,'frequency')
	#plt.arrow(70,10,0,-9,length_includes_head=True)
	#plt.arrow(0.765, 0.26, 0, -0.26, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	plt.arrow(0.768, 0.46, 0, -0.46, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	#plt.arrow(0.6857, 0.3, 0, -0.3, transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	#plt.text(1*10**(1),1500,'firing rate in')
	#plt.text(1*10**(1),500,'spiking state')
	plt.text(1*10**(1),3000,'firing rate in')
	plt.text(1*10**(1),1500,'spiking state')
	#plt.text(3*10**(-2),0.01,'firing rate')
	#plt.text(3*10**(-2),0.004,'in spiking state')
	#plt.text(8*10**(-7),2*d*timefac,'$2D_{eff}$')
	plt.text(10**(-7),2*d*timefac,'$2D_{eff}$')
	#plt.plot(10**(-5),2*d*timefac,'bx')
	#plt.arrow(0.10, 0.5785, -0.10, 0,transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
	plt.arrow(0.10, 0.294, -0.10, 0,transform=plot1.transAxes, length_includes_head=True,head_width=0.01,head_length=0.03)
#plt.plot(ax2,spectrum(ax2,1/(4.748*10**(-3)),13),label='theory')
#plt.plot(omega,background/T,'kx')
	#plt.legend(bbox_to_anchor=(0.4, 0.65))
	#plt.legend(loc='upper right')
	plt.legend()
#plt.plot(sax2,say2/T2,label='e6')
	#plt.savefig('specglue%s%sD=%.2fI=%.2f.pdf' %(mode,date,D*0.01,-0.1+z*0.02))
	#plot1.set_yticks(list(plot1.get_yticks()) +[2*d*timefac])
	plot1.spines['right'].set_visible(False)
	plot1.spines['top'].set_visible(False)
	plt.tight_layout()
	plt.savefig('specarinzelnew2.pdf')
	#plt.savefig('inapikanhopf%s2%sfourierD=%.2fI=%.2f.pdf' %(mode,date,D/100,43+z*0.25))	
