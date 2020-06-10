#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

from scipy.optimize import curve_fit

#matplotlib.rcParams.update({'font.size': 18})

def barrier(a,b,x):
	return a * x + b
def r(r0,U,D):
	return r0*np.exp(-U/D)
def snr(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,av,v0):
	return (qr(r0m,am,bm,cm,t,D)*(qr(r0p,ap,bp,cp,t,D)+qr(r0m,am,bm,cm,t,D))/(v0**2*qr(r0p,ap,bp,cp,t,D)))*((((2*ap*t+bp)-(2*am*t+bm))*qr(r0p,ap,bp,cp,t,D)*v0)/(D*(qr(r0p,ap,bp,cp,t,D)+qr(r0m,am,bm,cm,t,D)))+av)**2
def qbarrier(x,a,b,c):
	return a*x**2+b*x+c
def qr(r0,a,b,c,t,D):
	return r0*np.exp(-qbarrier(t,a,b,c)/D)
def deff(x,rp,rm,v0,av):
	return ((v0+av*x)**2*rp*rm)/((rp+rm)**3)
def fano(x,rp,rm,v0,av):
	return (2*(v0+av*x)*rp)/((rp+rm)**2)
def func(x, a, b):
	return a * np.exp(-b * x)
def deffana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return ((barrier(av,v0,t))**2*r(r0p,ap,bp,t,D)*r(r0m,am,bm,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**3)
def deffqana(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,av,v0):
	return ((barrier(av,v0,t))**2*qr(r0p,ap,bp,cp,t,D)*qr(r0m,am,bm,cm,t,D))/((qr(r0p,ap,bp,cp,t,D)+qr(r0m,am,bm,cm,t,D))**3)
def fanoana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return (2*barrier(av,v0,t)*r(r0p,ap,bp,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**2)
def vana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return barrier(av,v0,t)*r(r0m,am,bm,t,D)/(r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))
def vqana(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,av,v0):
	return barrier(av,v0,t)*qr(r0m,am,bm,cm,t,D)/(qr(r0p,ap,bp,cp,t,D)+qr(r0m,am,bm,cm,t,D))
def comp(x,b,c,d,e):
	return b*x**3+c*x**2+d*x+e
def comps(x,b,c,d):
	return 3*b*x**2+2*c*x+d
def fit(x,a,b):
	return a/x+b
def rprime(r0,r0s,Us,D):
	return (r0s/r0-Us/D)
def snrcor(r0p,r0m,up,um,ups,ums,D,av,v0,r0ps,r0ms):
	return ((av*r(r0m,um,D)*(r(r0m,um,D)+r(r0p,up,D))+v0*(rprime(r0m,r0ms,ums,D)-rprime(r0p,r0ps,ups,D))*r(r0p,up,D)*r(r0m,um,D))**2)/((r(r0p,up,D)+r(r0m,um,D))*v0**2*r(r0p,up,D)*r(r0m,um,D))

#date3='realrinzel25o'
#date2='realrinzel15ninv0'

date1='realrinzelrange26d1'
date2='realrinzelrange26d1'

istart=1
ivalues=10
#l=3
#D1=[]
#D3=[500]
#D2=[200,300]
#Dvar=[]
#D=D1+D2+D3+Dvar
#Da=np.array(D)
#btoeq=np.zeros((l,ivalues))
#eqtob=np.zeros((l,ivalues))
#params=np.zeros((4,ivalues))
#paramsav=np.zeros((4,ivalues))
#params2=np.zeros(4)
#paramsq=np.zeros(6)
#paramsqrate=np.zeros(6)
#for k2 in range(0,ivalues):
#	x=[]
#	ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/rate%s%d.txt' %(date3+date2,k2),'r')
#	for k4 in ratefile:
#		row=k4.split()
#		x.append(float(row[0]))
#	ax=np.array(x)
#	for k in range(0,l):
#		btoeq[k][k2]=1/ax[k]
#		eqtob[k][k2]=1/ax[k+l]
#

#xs=np.zeros(l)
#for b in range(0,l):
#	xs[b]=10/Da[b]
#for k2 in range(0,ivalues):
#	popt,pcov = curve_fit(func, xs, btoeq[:,k2])
#	params[0][k2]=popt[0]
#	params[1][k2]=popt[1]
#	popt,pcov = curve_fit(func, xs, eqtob[:,k2])
#	params[2][k2]=popt[0]
#	params[3][k2]=popt[1]
#rbte=np.mean(params[0,:])
#retb=np.mean(params[2,:])

#istart=1	
#xnew=np.arange(-22.5+istart,-22.5+istart+ivalues)*0.8

#barrier fit
#popt,pcov = curve_fit(qbarrier, xnew, params[1,:])
#paramsq[0]=popt[0]
#paramsq[1]=popt[1]
#paramsq[2]=popt[2]
#popt,pcov = curve_fit(qbarrier, xnew, params[3,:])
#paramsq[3]=popt[0]
#paramsq[4]=popt[1]
#paramsq[5]=popt[2]

#prefactor fit
#popt,pcov = curve_fit(qbarrier, xnew, params[0,:])
#paramsqrate[0]=popt[0]
#paramsqrate[1]=popt[1]
#paramsqrate[2]=popt[2]
#popt,pcov = curve_fit(qbarrier, xnew, params[2,:])
#paramsqrate[3]=popt[0]
#paramsqrate[4]=popt[1]
#paramsqrate[5]=popt[2]

#t1=np.arange(-18,-9,0.1)
#plt.figure()
#plt.plot(xnew,params[1,:],'go',label='burst to eq')
#plt.plot(xnew,params[3,:],'ro',label='eq to burst')
#plt.plot(t1,qbarrier(t1,paramsq[0],paramsq[1],paramsq[2]),'g')
#plt.plot(t1,qbarrier(t1,paramsq[3],paramsq[4],paramsq[5]),'r')
#plt.legend()
#plt.savefig('barrierinzelfit4.pdf')

#plt.figure()
#plt.plot(xnew,params[0,:],'go',label='burst to eq')
#plt.plot(xnew,params[2,:],'ro',label='eq to burst')
#plt.plot(t1,qbarrier(t1,paramsqrate[0],paramsqrate[1],paramsqrate[2]),'g')
#plt.plot(t1,qbarrier(t1,paramsqrate[3],paramsqrate[4],paramsqrate[5]),'r')
#plt.legend()
#plt.savefig('raterinzelfit4.pdf')

file=open('/home/richard/mastergit/pythonplots/arrhenius_analytics/detmocountparamparam.txt',"r")
col,colx=[],[]
for k in file:
	row=k.split()
	colx.append(float(row[0]))
	col.append(float(row[1]))
cola=np.array(col)
colxa=np.array(colx)

rbtoeq=np.zeros(ivalues)
ubtoeq=np.zeros(ivalues)
reqtob=np.zeros(ivalues)
ueqtob=np.zeros(ivalues)

for k2 in range(0,ivalues):
    x=[]
    ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/param%s%d.txt' %(date1+date2,k2),'r')
    for k4 in ratefile:
        row=k4.split()
        x.append(float(row[0]))
    ax=np.array(x)
    rbtoeq[k2]=x[0]
    ubtoeq[k2]=x[1]
    reqtob[k2]=x[2]
    ueqtob[k2]=x[3]

ups=np.zeros(ivalues-2)
ums=np.zeros(ivalues-2)
r0ms=np.zeros(ivalues-2)
r0ps=np.zeros(ivalues-2)
	
istep=0.8
xnew=np.arange(-21.25+istart,-21.25+istart+ivalues)*istep

for k3 in range(0,ivalues-2):
	ups[k3]=(ubtoeq[k3+2]-ubtoeq[k3])/(2*istep)
	ums[k3]=(ueqtob[k3+2]-ueqtob[k3])/(2*istep)
	r0ps[k3]=(rbtoeq[k3+2]-rbtoeq[k3])/(2*istep)
	r0ms[k3]=(reqtob[k3+2]-reqtob[k3])/(2*istep)

date='realrinzelrangelong26d1'
date1='realrinzelrange26d1'
date0='realrinzelrangeshort26d1'
date2='realrinzel13sig1'

D=[200]
D1=[250,300]
D0=[400,500]
D2=[]
Dto=D+D1+D0+D2
Dtot=np.array(Dto)
l2=len(Dtot)

points=1000000
length=500000
ivalues=10
istart=1
SNR=np.zeros((l2,ivalues))
scale=np.zeros((l2,ivalues))
xvec=np.zeros((l2,ivalues))
offset=np.zeros(l2,dtype=int)
ii=0

for c in D:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date,c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		if len(x) == 0:
			offset[ii]=offset[ii]+1
			continue
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
		scale[ii][z-istart-offset[ii]]=epsilon**2*T
		omegaind=round(omega*T)	
		#if c==200 and z==9:
		#	omegaind=2105	
		#omegaind=141	
		SNR[ii][z-istart-offset[ii]]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
		iv=open('/home/richard/outhome/d%s%d%d.txt' %(date,c,z),"r")
		for k in iv:
			row=k.split()
		xvec[ii][z-istart-offset[ii]]=float(row[0])
	ii=ii+1

for c1 in D1:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date1,c1,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		if len(x) == 0:
			offset[ii]=offset[ii]+1
			continue
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
		scale[ii][z-istart-offset[ii]]=epsilon**2*T
		omegaind=round(omega*T)		
		#omegaind=141
		SNR[ii][z-istart-offset[ii]]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
		iv=open('/home/richard/outhome/d%s%d%d.txt' %(date1,c1,z),"r")
		for k in iv:
			row=k.split()
		xvec[ii][z-istart-offset[ii]]=float(row[0])
	ii=ii+1

for c0 in D0:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date0,c0,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		if len(x) == 0:
			offset[ii]=offset[ii]+1
			continue
		ax=np.array(x)
		ay=np.array(y)
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date0,c0,z),"r")
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
		scale[ii][z-istart-offset[ii]]=epsilon**2*T
		omegaind=round(omega*T)		
		#omegaind=141
		SNR[ii][z-istart-offset[ii]]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
		iv=open('/home/richard/outhome/d%s%d%d.txt' %(date0,c0,z),"r")
		for k in iv:
			row=k.split()
		xvec[ii][z-istart-offset[ii]]=float(row[0])
	ii=ii+1

for c2 in D2:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date2,c2,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		if len(x) == 0:
			offset[ii]=offset[ii]+1
			continue
		ax=np.array(x)
		ay=np.array(y)
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
		scale[ii][z-istart-offset[ii]]=epsilon**2*T
		omegaind=round(omega*T)		
		#omegaind=141
		SNR[ii][z-istart-offset[ii]]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
		iv=open('/home/richard/outhome/d%s%d%d.txt' %(date2,c2,z),"r")
		for k in iv:
			row=k.split()
		xvec[ii][z-istart-offset[ii]]=float(row[0])
	ii=ii+1

#params=4
#dvalues=6
#b=np.zeros(dvalues)
#c=np.zeros(dvalues)
#d=np.zeros(dvalues)
#e=np.zeros(dvalues)
#ii=0
#xnew=np.arange(-0.19,0.31,0.01)
#eqfile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/detmocountparam5.txt','r')
#for k in eqfile:
#	col=[]
#	row=k.split()
#	for l in range(0,params):	
#		col.append(float(row[l]))
#	cola=np.array(col)
#	b[ii]=cola[0]
#	c[ii]=cola[1]
#	d[ii]=cola[2]
#	e[ii]=cola[3]
#	ii=ii+1
#eqfile.close() 
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
plt.ylabel('Signal-to-noise ratio SNR')

t=np.arange(-18,-8,0.1)
#xs=np.arange(-21.25+istart,-21.25+istart+ivalues)*0.8
xs=np.arange(-20+istart,-20+istart+ivalues)*0.6
plt.yscale('log')
#plt.xscale('log')
#plt.ylim(7*10**(-8),5*10**(-3))
plt.xlim(xnew[1]-0.4,xnew[ivalues-2]+0.4)
#plt.xlim(4*10**(-4),100)
#colorv=['r','y','c','g','k','b'] #6 colors
colorv=[ '#1f77b4', '#ff7f0e', '#2ca02c','#d62728','#9467bd'] 
#colorv=['y','g','b'] # 3 colors
#colorv=['y','c','g','k','b'] # 5 colors
#for n in range(0,1):
#	nl=round(ivalues-offset[n])
	#plt.plot(xvec[n,0:nl],(SNR[n,0:nl]-1)/scale[n,0:nl],colorv[n]+'o',label='D=%.2f' %(Dtot[n]*0.1))
#	plt.plot(xvec[n,0:nl],abs((SNR[n,0:nl]-1))/scale[n,0:nl],'o',label='D=%.0f*' %(Dtot[n]*0.1))
for n in range(0,l2):
	nl=round(ivalues-offset[n])
#	plt.plot(xvec[n,0:nl],(SNR[n,0:nl]-1)/scale[n,0:nl],colorv[n]+'o',label='D=%.2f' %(Dtot[n]*0.1))
	plt.plot(xvec[n,0:nl],abs((SNR[n,0:nl]-1))/scale[n,0:nl],label='D=%.0f' %(Dtot[n]*0.1))
#Dtot=[30,40,50,60]
#l2=len(Dtot)
#for n in range(0,l2):
	#plt.plot(xnew[1:ivalues-1],snrcor(rbtoeq[1:ivalues-1],reqtob[1:ivalues-1],ubtoeq[1:ivalues-1],ueqtob[1:ivalues-1],ups,ums,Dtot[n]*0.1,comps(xnew[1:ivalues-1],fit(Dtot[n]/10,colxa[0],cola[0]),fit(Dtot[n]/10,colxa[1],cola[1]),fit(Dtot[n]/10,colxa[2],cola[2])),comp(xnew[1:ivalues-1],fit(Dtot[n]/10,colxa[0],cola[0]),fit(Dtot[n]/10,colxa[1],cola[1]),fit(Dtot[n]/10,colxa[2],cola[2]),fit(Dtot[n]/10,colxa[3],cola[3])),r0ps,r0ms)/8,colorv[n])#,label='D=%.0f' %(Dtot[n]*0.1))
	#bv=b[2*n+1] # 3 plots
	#cv=c[2*n+1]
	#dv=d[2*n+1]
	#ev=e[2*n+1]
	#bv=b[n+1] # 5
	#cv=c[n+1]
	#dv=d[n+1]
	#ev=e[n+1]	#plt.plot(t,snr(rbte,retb,paramsq[0],paramsq[3],paramsq[1],paramsq[4],paramsq[2],paramsq[5],t,Dtot[n]*0.1,comps(t,b[n+1],c[n+1],d[n+1]),comp(t,b[n+1],c[n+1],d[n+1],e[n+1]))/8,colorv[n])	
	#plt.plot(t,snr(rbte,retb,paramsq[0],paramsq[3],paramsq[1],paramsq[4],paramsq[2],paramsq[5],t,Dtot[n]*0.1,comps(t,bv,cv,dv),comp(t,bv,cv,dv,ev))/8,colorv[n])
	#plt.plot(t,snr(qbarrier(t,paramsqrate[0],paramsqrate[1],paramsqrate[2]),qbarrier(t,paramsqrate[3],paramsqrate[4],paramsqrate[5]),paramsq[0],paramsq[3],paramsq[1],paramsq[4],paramsq[2],paramsq[5],t,Dtot[n]*0.1,comps(t,bv,cv,dv),comp(t,bv,cv,dv,ev))/8,colorv[n])
#plt.plot(t,snr(qbarrier(t,paramsqrate[0],paramsqrate[1],paramsqrate[2]),qbarrier(t,paramsqrate[3],paramsqrate[4],paramsqrate[5]),paramsq[0],paramsq[3],paramsq[1],paramsq[4],paramsq[2],paramsq[5],t,Dtot[3]*0.1,comps(t,b[5],c[5],d[5]),comp(t,b[5],c[5],d[5],e[5]))/8,colorv[3])
#plt.plot([0.163, 0.163], [10**(-7), 10], color='black', linestyle='-')
#plt.plot([-0.02, -0.02], [10**(-7), 10], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [0,2,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
plt.plot([-10.8, -10.8], [10**(-8), 10**(-3)], color='black', linestyle='-',label='$I_{crit}$')
plt.plot([-13, -13], [10**(-8), 10**(-3)], color='black', linestyle='--',label='$I_{max}$')
#plt.plot([-10.8, -10.8], [10**(-38), 10**(5)], color='black', linestyle='-',label='$I_{crit}$')
plt.legend()
plt.tight_layout()
#plt.savefig('snrinzelrange26dcompletecritnofit3big.pdf')
#plt.savefig('snrinzelpred2big.pdf')
plt.savefig('snrinzelonlycritmax.pdf')
