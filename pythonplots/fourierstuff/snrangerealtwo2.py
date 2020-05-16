#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

from scipy.optimize import curve_fit

#matplotlib.rcParams.update({'font.size': 14})

def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
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


date3='realfast11jjem2sh'
date2='realfast19jjem2st'

ivalues=20
l=5
D1=[35]
D3=[40,50]
D2=[45]
Dvar=[30]
D=D1+D2+D3+Dvar
Da=np.array(D)
btoeq=np.zeros((l,ivalues))
eqtob=np.zeros((l,ivalues))
params=np.zeros((4,ivalues))
paramsav=np.zeros((4,ivalues))
params2=np.zeros(4)
paramsq=np.zeros(6)
for k2 in range(0,ivalues):
	x=[]
	ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/rate%s%d.txt' %('new'+date3+'new'+date2,k2),'r')
	for k4 in ratefile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	for k in range(0,l):
		btoeq[k][k2]=1/ax[k]
		eqtob[k][k2]=1/ax[k+l]

av=0.013
v0=0.0637
xs=np.zeros(l)
for b in range(0,l):
	xs[b]=100/Da[b]
for k2 in range(0,ivalues):
	popt,pcov = curve_fit(func, xs, btoeq[:,k2])
	params[0][k2]=popt[0]
	params[1][k2]=popt[1]
	popt,pcov = curve_fit(func, xs, eqtob[:,k2])
	params[2][k2]=popt[0]
	params[3][k2]=popt[1]
rbte=np.mean(params[0,:])
retb=np.mean(params[2,:])

istart=1	
xnew=np.arange(-5+istart,-5+istart+ivalues)*0.02

popt,pcov = curve_fit(qbarrier, xnew, params[1,:])
paramsq[0]=popt[0]
paramsq[1]=popt[1]
paramsq[2]=popt[2]
popt,pcov = curve_fit(qbarrier, xnew, params[3,:])
paramsq[3]=popt[0]
paramsq[4]=popt[1]
paramsq[5]=popt[2]


bp = 1.74
ap = 5.64
r0p = 0.0075
bm = 3.15
am = -10.76
r0m = 0.012

date=['realfast13aem2n4','realfast13aem2n4','realfast9aem2sh','realfast9aem2sh']

D=[25,30,35,45]
l=len(D)

points=1000000
length=500000
ivalues=20
SNR=np.zeros((l,ivalues))
scale=np.zeros((l,ivalues))
ii=0

for m in range(0,l):
	c=D[m]
	for z in range(1,21):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date[m],c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date[m],c,z),"r")
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
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Signal-to-noise ratio SNR')

t=np.arange(-0.1,0.3,0.01)
xs=np.arange(-0.08,0.32,0.02)
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
#colorv=['g','y','b','r','c']
for n in range(0,l):
	plt.plot(xs,abs(SNR[n,:]-1)/scale[n,:],label='D=%.2f' %(D[n]*0.01))
#for n in range(0,l):	
#	plt.plot(t,snr(rbte,retb,paramsq[0],paramsq[3],paramsq[1],paramsq[4],paramsq[2],paramsq[5],t,Dtot[n]*0.01,av,v0)/8,colorv[n])
plt.plot([0.163, 0.163], [5*10**(-5), 5*10**(-2)], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot([-0.02, -0.02], [10**(-7), 10], color='black', linestyle='-')
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [0,2,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
plt.legend()
plt.tight_layout()
plt.savefig('snrealonly2crit4.pdf')
