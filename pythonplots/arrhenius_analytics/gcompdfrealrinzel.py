#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

from scipy.optimize import curve_fit

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
def deffqana(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,v0):
	return (v0**2*qr(r0p,ap,bp,cp,t,D)*qr(r0m,am,bm,cm,t,D))/((qr(r0p,ap,bp,cp,t,D)+qr(r0m,am,bm,cm,t,D))**3)
def fanoana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return (2*barrier(av,v0,t)*r(r0p,ap,bp,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**2)
def fanoqana(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,v0):
	return 2*deffqana(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,v0)/vqana(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,v0)
def vana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return barrier(av,v0,t)*r(r0m,am,bm,t,D)/(r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))
def vqana(r0p,r0m,ap,am,bp,bm,cp,cm,t,D,v0):
	return v0*qr(r0m,am,bm,cm,t,D)/(qr(r0p,ap,bp,cp,t,D)+qr(r0m,am,bm,cm,t,D))
def comp(x,b,c,d,e):
	return b*x**3+c*x**2+d*x+e
def comps(x,b,c,d):
	return 3*b*x**2+2*c*x+d

date3='realrinzel25o'
date2='realrinzel15ninv0'

ivalues=11
l=3
D1=[]
D3=[200,300,500]
D2=[]
Dvar=[]
D=D1+D2+D3+Dvar
Da=np.array(D)
btoeq=np.zeros((l,ivalues))
eqtob=np.zeros((l,ivalues))
params=np.zeros((4,ivalues))
paramsav=np.zeros((4,ivalues))
params2=np.zeros(4)
paramsq=np.zeros(6)
paramsqrate=np.zeros(6)

for k2 in range(0,ivalues):
	x=[]
	ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/rate%s%d.txt' %(date3+date2,k2),'r')
	for k4 in ratefile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	for k in range(0,l):
		btoeq[k][k2]=1/ax[k]
		eqtob[k][k2]=1/ax[k+l]

av=0.1/12
v0=0.47
xs=np.zeros(l)
for b in range(0,l):
	xs[b]=10/Da[b]
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
xnew=np.arange(-22.5+istart,-22.5+istart+ivalues)*0.8

#barrier fit
popt,pcov = curve_fit(qbarrier, xnew, params[1,:])
paramsq[0]=popt[0]
paramsq[1]=popt[1]
paramsq[2]=popt[2]
popt,pcov = curve_fit(qbarrier, xnew, params[3,:])
paramsq[3]=popt[0]
paramsq[4]=popt[1]
paramsq[5]=popt[2]

#prefactor fit
popt,pcov = curve_fit(qbarrier, xnew, params[0,:])
paramsqrate[0]=popt[0]
paramsqrate[1]=popt[1]
paramsqrate[2]=popt[2]
popt,pcov = curve_fit(qbarrier, xnew, params[2,:])
paramsqrate[3]=popt[0]
paramsqrate[4]=popt[1]
paramsqrate[5]=popt[2]

vec=np.zeros((l,ivalues))
vecx=np.zeros((l,ivalues))
ii=0

date='realrinzel25o'
date1='realrinzel15ninv0'
D=[200,300,500]
D1=[]
Dtot=D+D1
l=len(D)+len(D1)

for x in D:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D1:
	col,colx=[],[]	
	for y in range(2,10):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,8):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1


dvalues=6
dparams=4
b=np.zeros(dvalues)
c=np.zeros(dvalues)
d=np.zeros(dvalues)
e=np.zeros(dvalues)
ii=0
#xnew=np.arange(-0.19,0.31,0.01)
eqfile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/detmocountparam5.txt','r')
for k in eqfile:
	col=[]
	row=k.split()
	for l in range(0,dparams):	
		col.append(float(row[l]))
	cola=np.array(col)
	b[ii]=cola[0]
	c[ii]=cola[1]
	d[ii]=cola[2]
	e[ii]=cola[3]
	ii=ii+1
eqfile.close() 
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
plt.ylabel('firing rate')

t=np.arange(-18,-9,0.1)
xs=np.arange(-22.5+istart,-22.5+istart+ivalues)*0.8
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
colorv=['y','g','b','r','c']

for n in range(0,l):
	plt.plot(vecx[n,:],vec[n,:],colorv[n]+'o',label='D=%.2f' %(Dtot[n]*0.1))
for n in range(0,l):
	bv=b[2*n+1]
	cv=c[2*n+1]
	dv=d[2*n+1]
	ev=e[2*n+1]	#plt.plot(t,snr(rbte,retb,paramsq[0],paramsq[3],paramsq[1],paramsq[4],paramsq[2],paramsq[5],t,Dtot[n]*0.1,comps(t,b[n+1],c[n+1],d[n+1]),comp(t,b[n+1],c[n+1],d[n+1],e[n+1]))/8,colorv[n])	
	plt.plot(t,deffqana(qbarrier(t,paramsqrate[0],paramsqrate[1],paramsqrate[2]),qbarrier(t,paramsqrate[3],paramsqrate[4],paramsqrate[5]),paramsq[0],paramsq[3],paramsq[1],paramsq[4],paramsq[2],paramsq[5],t,Dtot[n]*0.1,comp(t,bv,cv,dv,ev)),colorv[n])
#plt.plot([0.163, 0.163], [10**(-7), 10], color='black', linestyle='-')
#plt.plot([-0.02, -0.02], [10**(-7), 10], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [0,2,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
plt.legend()
plt.savefig('dcomprate25o6j.pdf')
