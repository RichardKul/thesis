#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def barrier(a,b,x):
	return a * x + b
def qbarrier(x,a,b,c):
	return a*x**2+b*x+c
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def deff(x,rp,rm,v0,av):
	return ((v0+av*x)**2*rp*rm)/((rp+rm)**3)
def func(x, a, b):
	return a * np.exp(-b * x)
def deffana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return ((barrier(av,v0,t))**2*r(r0p,ap,bp,t,D)*r(r0m,am,bm,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**3)
def fanoana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return (2*barrier(av,v0,t)*r(r0p,ap,bp,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**2)
def vana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return barrier(av,v0,t)*r(r0m,am,bm,t,D)/(r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))


D1=[35]
D3=[40,50]
D2=[45]
Dvar=[30]
D=D1+D2+D3+Dvar
Da=np.array(D)
l1=len(D1)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l1+l2+l3+lvar
date1='realfast11jjem2sh'
date3='realfast11jjem2st'
date2='realfast19jjem2st'
datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
yvar=[4,13,3]
yvalues=len(yvar)

istart=1
ivalues=20
epsilon = 0.00001

btoeq=np.zeros((l,ivalues))
eqtob=np.zeros((l,ivalues))
params=np.zeros((4,ivalues))
paramsav=np.zeros((4,ivalues))
params2=np.zeros(4)
paramsq=np.zeros(6)
for k2 in range(0,ivalues):
	x=[]
	ratefile = open('rate%s%d.txt' %('new'+date1+'new'+date2,k2),'r')
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


for k2 in range(0,ivalues):
	popt,pcov = curve_fit(func, xs, btoeq[:,k2],bounds=((rbte-epsilon,-np.inf), (rbte+epsilon,np.inf)))
	paramsav[0][k2]=popt[0]
	paramsav[1][k2]=popt[1]
	popt,pcov = curve_fit(func, xs, eqtob[:,k2],bounds=((retb-epsilon,-np.inf), (retb+epsilon,np.inf)))
	paramsav[2][k2]=popt[0]
	paramsav[3][k2]=popt[1]
	
xnew=np.arange(-5+istart,-5+istart+ivalues)*0.02

popt,pcov = curve_fit(qbarrier, xnew, params[1,:])
paramsq[0]=popt[0]
paramsq[1]=popt[1]
paramsq[2]=popt[2]
popt,pcov = curve_fit(qbarrier, xnew, params[3,:])
paramsq[3]=popt[0]
paramsq[4]=popt[1]
paramsq[5]=popt[2]

popt,pcov = curve_fit(barrier, xnew, paramsav[1,:])
params2[0]=popt[0]
params2[1]=popt[1]
popt,pcov = curve_fit(barrier, xnew, paramsav[3,:])
params2[2]=popt[0]
params2[3]=popt[1]
eqfile2 = open('paramsdfnew%s.txt'%(date1+date2),'w')
for k4 in range(0,2): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.write('%.6f\n'%rbte) 
for k4 in range(2,4): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.write('%.6f\n'%retb) 
eqfile2.close() 
	
paramfile = open('parambarrier%s.txt' %('new'+date1+'new'+date2),'r')
xx=[]
for k4 in paramfile:
	row=k4.split()
	xx.append(float(row[0]))
axx=np.array(xx)



ii=0

istart1=1
ivalues1=20

vec=np.zeros((l,ivalues1))

for x in D1:
	col1,colx1=[],[]	
	for y in range(istart1,istart1+ivalues1):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues1):
		vec[ii][z]=cola1[z]
	ii=ii+1

istart2=1
ivalues2=20

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

xs=np.arange(-5+istart1,-5+istart1+ivalues1)*0.02
pv=np.arange(0,ivalues)

plt.xlabel('bias current I')
plt.ylabel('$D_{eff}$')
#plt.xlim(0,3.5)
#plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')

ap=ax[0]
am=ax[3]
bp=ax[1]
bm=ax[4]
r0p=rbte
r0m=retb

colorv=['y','g','b','r','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#for n in range(0,l):
#	plt.plot(t,deffana(axx[2],axx[5],axx[0],axx[3],axx[1],axx[4],t,Da[n]/100,av,v0),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot(xs,vec[n,:],colorv[n]+'o',label='D=%.2f' %(Da[n]/100))
for n in range(0,l):
	plt.plot((pv-4)*0.02,deff((pv-4)*0.02,btoeq[n][pv],eqtob[n][pv],v0,av),colorv[n])
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')


handles, labels = plt.gca().get_legend_handles_labels()
order = [4,0,2,1,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('dcompdfpwnew%s.pdf' %(date1+date2))

g=np.zeros((l,ivalues1))
ii=0

for x in D1:
	col1,colx1=[],[]	
	for y in range(istart1,istart1+ivalues1):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues1):
		g[ii][z]=cola1[z]
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
		g[ii][z]=cola2[z]
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
		g[ii][z]=cola3[z]
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
		g[ii][z]=colavar[z]
	ii=ii+1

plt.figure()
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.yscale('log')

colorv=['y','g','b','r','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,vana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot(t,vana(axx[2],axx[5],axx[0],axx[3],axx[1],axx[4],t,Da[n]/100,av,v0),colorv[n])
for n in range(0,l):
	plt.plot(xs,g[n,:],colorv[n]+'o',label='D=%.2f' %(Da[n]/100))
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

handles, labels = plt.gca().get_legend_handles_labels()
order = [4,0,2,1,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('gcompdfnew%s.pdf' %(date1+date2))


fano=np.zeros((l,ivalues1))
ii=0

for x in D1:
	col1,colx1=[],[]	
	for y in range(istart1,istart1+ivalues1):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues1):
		fano[ii][z]=cola1[z]
	ii=ii+1

for x in D2:
	col2,colx2=[],[]	
	for y in range(istart2,istart2+ivalues2):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx2.append(float(row[0]))
			col2.append(float(row[1]))
	colxa2=np.array(colx2)
	cola2=np.array(col2)
	for z in range(0,ivalues2):
		fano[ii][z]=cola2[z]
	ii=ii+1

for x in D3:
	col3,colx3=[],[]	
	for y in range(istart3,istart3+ivalues3):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx3.append(float(row[0]))
			col3.append(float(row[1]))
	colxa3=np.array(colx3)
	cola3=np.array(col3)
	for z in range(0,ivalues3):
		fano[ii][z]=cola3[z]
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
			file=open('/home/richard/outhome/f%s%d%d.txt' % (datevar[kk],x,y),"r")
			for k in file:
				row=k.split()
				colxvar.append(float(row[0]))
				colvar.append(float(row[1]))
	colxavar=np.array(colxvar)
	colavar=np.array(colvar)
	for z in range(0,20):
		fano[ii][z]=colavar[z]
	ii=ii+1

plt.figure()
plt.xlabel('bias current I')
plt.ylabel('Fano factor')
plt.yscale('log')

colorv=['y','g','b','r','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
for n in range(0,l):
	plt.plot(t,fanoana(axx[2],axx[5],axx[0],axx[3],axx[1],axx[4],t,Da[n]/100,av,v0),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot(xs,fano[n,:],colorv[n]+'o',label='D=%.2f' %(Da[n]/100))
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')


handles, labels = plt.gca().get_legend_handles_labels()
order = [4,0,2,1,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('fcompdfnew%s.pdf' %(date1+date2))

t=np.arange(-0.1,0.3,0.001)
plt.figure()
plt.xlabel('bias current')
plt.ylabel('potential barrier')
plt.plot(t,barrier(t,params2[0],params2[1]),'y')
plt.plot(t,barrier(t,params2[2],params2[3]),'y')
plt.plot(t,2*barrier(t,params2[0],params2[1]),'y')
plt.plot(t,2*barrier(t,params2[2],params2[3]),'y')
plt.plot(t,qbarrier(t,paramsq[0],paramsq[1],paramsq[2]),'y')
plt.plot(t,qbarrier(t,paramsq[3],paramsq[4],paramsq[5]),'y')
plt.plot(t,2*qbarrier(t,paramsq[0],paramsq[1],paramsq[2]),'y')
plt.plot(t,2*qbarrier(t,paramsq[3],paramsq[4],paramsq[5]),'y')
plt.plot(xs,params[1,:],label='burst to eq')
plt.plot(xs,params[3,:],label='eq to burst')
plt.plot(xs,2*params[1,:],label='2x burst to eq')
plt.plot(xs,2*params[3,:],label='2x eq to burst')
plt.legend()
plt.savefig('barriercompdfnew%s.pdf' %(date1+date2))
