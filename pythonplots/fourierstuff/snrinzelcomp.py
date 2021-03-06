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
def drdi(t,r0,a,b):
	return (a*r0)/((1+np.exp(-barrier(a,b,t)))*(1+np.exp(barrier(a,b,t))))
def snr(r0p,r0m,ap,am,bp,bm,t,D):
	return (ap-am)**2/(D**2*(r(r0p,ap,bp,t,D)**(-1)+r(r0m,am,bm,t,D)**(-1)))
def snrtheo(t,r0,a,b,f):
	return drdi(t,r0,a,b)**2/(8*f**2)
def snrmeas(drdr,f):
	return (drdr**2/(8*f**2)) 


f=1



date='realrinzel25o'
date1='realfast9aem2sh'
D=[200,300,500]
D1=[]
Dtot=D+D1
l=len(D)+len(D1)

points=1000000
length=500000
ivalues=10
SNR=np.zeros((l,ivalues))
scale=np.zeros((l,ivalues))
ii=0

#for c in D:
#	for z in range(1,21):
#		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date,c,z),"r")
#		x,y=[],[]
#		for k in file:
#			row=k.split()
#			x.append(float(row[0]))
#			y.append(float(row[1]))
#		ax=np.array(x)
#		ay=np.array(y)
#		param=open('/home/richard/outhome/param%s%d%d.txt' %(date,c,z),"r")
#		ll=0
#		name,value=[],[]
#		for k in param:
#			row=k.split()
#			lp=len(row)
#			if ll<1:
#				for jj in range(0,lp):
#					name.append(row[jj])
#			else:
#				for kk in range(0,lp):
#					value.append(float(row[kk]))
#			ll=ll+1
#		dt=value[name.index('dt')]
#		N=value[name.index('N')]-value[name.index('Neq')]
#		repetitions=value[name.index('repetitions')]
#		epsilon=value[name.index('epsilon')]
#		omega=value[name.index('omega')]
#		T=N*repetitions*dt
#		scale[ii][z-1]=epsilon**2*T
#		omegaind=round(omega*T)		
#		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
#	ii=ii+1
#
#for c1 in D1:
#	for z in range(1,21):
#		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date1,c1,z),"r")
#		x,y=[],[]
#		for k in file:
#			row=k.split()
#			x.append(float(row[0]))
#			y.append(float(row[1]))
#		ax=np.array(x)
#		ay=np.array(y)
#		param=open('/home/richard/outhome/param%s%d%d.txt' %(date1,c1,z),"r")
#		ll=0
#		name,value=[],[]
#		for k in param:
#			row=k.split()
#			lp=len(row)
#			if ll<1:
#				for jj in range(0,lp):
#					name.append(row[jj])
#			else:
#				for kk in range(0,lp):
#					value.append(float(row[kk]))
#			ll=ll+1
#		dt=value[name.index('dt')]
#		N=value[name.index('N')]-value[name.index('Neq')]
#		repetitions=value[name.index('repetitions')]
#		epsilon=value[name.index('epsilon')]
#		omega=value[name.index('omega')]
#		T=N*repetitions*dt
#		scale[ii][z-1]=epsilon**2*T
#		omegaind=round(omega*T)		
#		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
#	ii=ii+1



D0=[200,300,500]
D3=[]
D2=[]
Dvar=[]
D=D0+Dvar+D3+D2
Da=np.array(D)
l0=len(D0)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l2+l3+lvar+l0
date3='realfast11jjem2sh'
date2='realfast19jjem2st'
datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
yvar=[4,13,3]
yvalues=len(yvar)


ii=0

istart0=1
ivalues0=10
istart2=1
ivalues2=10
vec=np.zeros((l,ivalues))

for x in D0:
	col0,colx0=[],[]	
	for y in range(istart0,istart0+ivalues0):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,ivalues0):
		vec[ii][z]=cola0[z]
	ii=ii+1


for x in Dvar:
	colvar,colxvar=[],[]
	for kk in range(0,yvalues):
		ystart=1
		jj=0
		while jj < kk:
			ystart=ystart+yvar[jj]
			jj=jj+1
		for y in range(ystart,ystart+yvar[kk]):
			file=open('/home/richard/outhome/d%s%d%d.txt' % (datevar[kk],x,y),"r")
			for k in file:
				row=k.split()
				colxvar.append(float(row[0]))
				colvar.append(float(row[1]))
	colxavar=np.array(colxvar)
	colavar=np.array(colvar)
	for z in range(0,20):
		vec[ii][z]=colavar[z]
	ii=ii+1

istart3=1
ivalues3=20

for x in D3:
	col3,colx3=[],[]	
	for y in range(istart3,istart3+ivalues3):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx3.append(float(row[0]))
			col3.append(float(row[1]))
	colxa3=np.array(colx3)
	cola3=np.array(col3)
	for z in range(0,ivalues3):
		vec[ii][z]=cola3[z]
	ii=ii+1

for x in D2:
	col2,colx2=[],[]	
	for y in range(istart2,istart2+ivalues2):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx2.append(float(row[0]))
			col2.append(float(row[1]))
	colxa2=np.array(colx2)
	cola2=np.array(col2)
	for z in range(0,ivalues2):
		vec[ii][z]=cola2[z]
	ii=ii+1

istep=0.8

ii=0
gvec=np.zeros((l,ivalues))
#drdr=np.zeros((l,ivalues-2))

D0,D2,D3,Dvar=[],[],[],[]
for x in D0:
	col0,colx0=[],[]	
	for y in range(istart0,istart0+ivalues0):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,ivalues0):
		gvec[ii][z]=cola0[z]
	for z in range(0,8):
		drdr[ii][z]=(cola0[z+2]-cola0[z])/(2*istep)
	ii=ii+1

for x in Dvar:
	colvar,colxvar=[],[]
	for kk in range(0,yvalues):
		ystart=1
		jj=0
		while jj < kk:
			ystart=ystart+yvar[jj]
			jj=jj+1
		for y in range(ystart,ystart+yvar[kk]):
			file=open('/home/richard/outhome/g%s%d%d.txt' % (datevar[kk],x,y),"r")
			for k in file:
				row=k.split()
				colxvar.append(float(row[0]))
				colvar.append(float(row[1]))
	colxavar=np.array(colxvar)
	colavar=np.array(colvar)
	for z in range(0,20):
		gvec[ii][z]=colavar[z]
	for z in range(0,18):
		drdr[ii][z]=(colavar[z+2]-colavar[z])/(2*istep)
	ii=ii+1


for x in D3:
	col3,colx3=[],[]	
	for y in range(istart3,istart3+ivalues3):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx3.append(float(row[0]))
			col3.append(float(row[1]))
	colxa3=np.array(colx3)
	cola3=np.array(col3)
	for z in range(0,ivalues3):
		gvec[ii][z]=cola3[z]
	for z in range(0,18):
		drdr[ii][z]=(cola3[z+2]-cola3[z])/(2*istep)
	ii=ii+1

for x in D2:
	col2,colx2=[],[]	
	for y in range(istart2,istart2+ivalues2):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx2.append(float(row[0]))
			col2.append(float(row[1]))
	colxa2=np.array(colx2)
	cola2=np.array(col2)
	for z in range(0,ivalues2):
		gvec[ii][z]=cola2[z]
	for z in range(0,18):
		drdr[ii][z]=(cola2[z+2]-cola2[z])/(2*istep)
	ii=ii+1

istep=0.3
drdr=np.zeros((l,10))
drdrlong=np.zeros((l,50))
d2z,d3z,d5z=[],[],[]
filerate=open('/home/richard/outhome/countrinzelrate.txt',"r")
for k in filerate:
	row=k.split()
	d2z.append(float(row[2]))
	d3z.append(float(row[4]))
	d5z.append(float(row[6]))
d2=np.array(d2z)/10000000
d3=np.array(d3z)/10000000
d5=np.array(d5z)/10000000
for z in range(0,10):
	drdr[0][z]=abs((d2[round(9.5+(8/3)*(z+1))+1]-d2[round(9.5+(8/3)*(z+1))])/istep)
	drdr[1][z]=abs((d3[round(9.5+(8/3)*(z+1))+1]-d3[round(9.5+(8/3)*(z+1))])/istep)
	drdr[2][z]=abs((d5[round(9.5+(8/3)*(z+1))+1]-d5[round(9.5+(8/3)*(z+1))])/istep)
for z in range(0,50):
	drdrlong[0][z]=abs((d2[z+1]-d2[z])/istep)
	drdrlong[1][z]=abs((d3[z+1]-d3[z])/istep)
	drdrlong[2][z]=abs((d5[z+1]-d5[z])/istep)


xsold=np.arange(-69.5+1,-69.5+51)*0.3
plt.figure()
plt.xlabel('bias current I')
plt.ylabel('derivative of firing rate')
plt.yscale('log')
colorv=['y','g','b','r','c']
for n in range(0,l):
	plt.plot(xsold,drdrlong[n,:],colorv[n],label='D=%.0f' %Da[n])
plt.legend()
plt.savefig('drdirinzellong.pdf')

date='realrinzel20ninv0'
date2='realrinzel20ninv1'
date1='realrinzel15ninv0'
date3='realrinzel15ninv1'
D=[200,300]
D1=[500]
D2=[]
D3=[]
Dtot=D+D1+D2+D3
l=len(D)+len(D1)+len(D2)+len(D3)

points=1000000
length=500000
istart=1
ivalues=10
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
xs=np.arange(-16.2,-9,0.8)
xsh=np.arange(-17.2,-9.2,0.8)
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
plt.ylim(10**(-8),10**(-3))
colorv=['y','g','b','r','c']
for n in range(0,l):	
	plt.plot(xsh,snrmeas(drdr[n,:],f)/vec[n,:],colorv[n],label='D=%.2f' %(Dtot[n]/10))
for n in range(0,l):
	plt.plot(xs,abs((SNR[n,0:9]-1)/scale[n,0:9]),colorv[n]+'o')
#plt.plot([0.163, 0.163], [2*10**(-5), 6*10**(-2)], color='black', linestyle='-')
#plt.plot([-0.02, -0.02], [2*10**(-5), 6*10**(-2)], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [0,2,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
plt.legend()
plt.savefig('snrangerealrinzelcomplong.pdf')
