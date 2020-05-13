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
def v(x,rp,rm,v0,av):
	return ((v0+av*x)*rm)/(rp+rm)
def func(x, a, b):
	return a * np.exp(-b * x)
def deffana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return ((barrier(av,v0,t))**2*r(r0p,ap,bp,t,D)*r(r0m,am,bm,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**3)
def fanoana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return (2*barrier(av,v0,t)*r(r0p,ap,bp,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**2)
def vana(r0p,r0m,ap,am,bp,bm,t,D,av,v0):
	return barrier(av,v0,t)*r(r0m,am,bm,t,D)/(r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))
def fanofun(x,rp,rm,v0,av):
	return 2*deff(x,rp,rm,v0,av)/v(x,rp,rm,v0,av)
def comp(x,b,c,d,e):
	return b*x**3+c*x**2+d*x+e
def fit(x,a,b):
	return a/x+b

timefac=1000 

D2=[250,300]
D1=[200]
D3=[400,500]
Dvar=[]
D=D1+D2+D3+Dvar
Da=np.array(D)
l1=len(D1)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l1+l2+l3+lvar
date1='realrinzelrange26d1'
date3=''
date2='realrinzelrange26d1'
datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
yvar=[4,13,3]
yvalues=len(yvar)

istart=1
ivalues=10
epsilon = 0.00001


#params=4
#dvalues=6
#b=np.zeros(dvalues)
#c=np.zeros(dvalues)
#d=np.zeros(dvalues)
#e=np.zeros(dvalues)
#ii=0
#eqfile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/detmocountparam5.txt','r')
#for k in eqfile:
#	col=[]
#	row=k.split()
#	for lll in range(0,params):	
#		col.append(float(row[lll]))
#	cola=np.array(col)
#	b[ii]=cola[0]
#	c[ii]=cola[1]
#	d[ii]=cola[2]
#	e[ii]=cola[3]
#	ii=ii+1
#eqfile.close() 

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

def rred(r0,U,D):
	return r0*np.exp(-U/D)
def deffred(rp,up,rm,um,v,D):
	return (v**2*rred(rp,up,D)*rred(rm,um,D))/((rred(rp,up,D)+rred(rm,um,D))**3)
def vred(rp,up,rm,um,v,D):
	return v*rred(rm,um,D)/(rred(rp,up,D)+rred(rm,um,D))
def fanored(rp,up,rm,um,v,D):
	return 2*deffred(rp,up,rm,um,v,D)/vred(rp,up,rm,um,v,D)

date2='realrinzelrange26d1'
date3='realrinzelrangeshort26d1'
date1='realrinzelrangelong26d1'
ii=0

istart1=1
ivalues1=10

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
ivalues2=10

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
ivalues3=10

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

xs=np.arange(-21.25+istart1,-21.25+istart1+ivalues1)*0.8
pv=np.arange(0,ivalues)

plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ [$10^3s^{-1}$]')
#plt.xlim(0,3.5)
#plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')

colorv=[ '#1f77b4', '#ff7f0e', '#2ca02c','#d62728','#9467bd']
#for n in range(0,l):
#	plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#for n in range(0,l):
#	plt.plot(t,deffana(axx[2],axx[5],axx[0],axx[3],axx[1],axx[4],t,Da[n]/100,av,v0),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,1):
		plt.plot(xs,vec[n,:]*timefac,'o',color=colorv[n],label='D=%.0f*' %(Da[n]/10))
for n in range(1,l):
		plt.plot(xs,vec[n,:]*timefac,'o',color=colorv[n],label='D=%.0f' %(Da[n]/10))
for n in range(0,l):
	#bv=b[n+1]
	#cv=c[n+1]
	#dv=d[n+1]
	#ev=e[n+1]
	plt.plot((pv-20.25)*0.8,deffred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],comp((pv-20.25)*0.8,fit(Da[n]/10,colxa[0],cola[0]),fit(Da[n]/10,colxa[1],cola[1]),fit(Da[n]/10,colxa[2],cola[2]),fit(Da[n]/10,colxa[3],cola[3])),Da[n]/10)*timefac,colorv[n])
	
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')


#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('dcompdfpwnew2%s.pdf' %(date1+date2))

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
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('firing rate r [$10^3s^{-1}$]')
#plt.yscale('log')

colorv=[ '#1f77b4', '#ff7f0e', '#2ca02c','#d62728','#9467bd']
#for n in range(0,l):
#	plt.plot(t,vana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,1):
		plt.plot(xs,g[n,:]*timefac,'o',color=colorv[n],label='D=%.0f*' %(Da[n]/10))
for n in range(1,l):
		plt.plot(xs,g[n,:]*timefac,'o',color=colorv[n],label='D=%.0f' %(Da[n]/10))
for n in range(0,l):
	#bv=b[n+1]
	#cv=c[n+1]
	#dv=d[n+1]
	#ev=e[n+1]
	plt.plot((pv-20.25)*0.8,vred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],comp((pv-20.25)*0.8,fit(Da[n]/10,colxa[0],cola[0]),fit(Da[n]/10,colxa[1],cola[1]),fit(Da[n]/10,colxa[2],cola[2]),fit(Da[n]/10,colxa[3],cola[3])),Da[n]/10)*timefac,colorv[n])

#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('gcompdfpwnew2%s.pdf' %(date1+date2))


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
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor')
plt.yscale('log')

colorv=[ '#1f77b4', '#ff7f0e', '#2ca02c','#d62728','#9467bd']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,1):
		plt.plot(xs,fano[n,:],'o',color=colorv[n],label='D=%.0f*' %(Da[n]/10))
for n in range(1,l):
		plt.plot(xs,fano[n,:],'o',color=colorv[n],label='D=%.0f' %(Da[n]/10))
for n in range(0,l):
	#bv=b[n+1]
	#cv=c[n+1]
	#dv=d[n+1]
	#ev=e[n+1]
	plt.plot((pv-20.25)*0.8,fanored(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],comp((pv-20.25)*0.8,fit(Da[n]/10,colxa[0],cola[0]),fit(Da[n]/10,colxa[1],cola[1]),fit(Da[n]/10,colxa[2],cola[2]),fit(Da[n]/10,colxa[3],cola[3])),Da[n]/10),colorv[n])
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')


#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('fcompdfpwnew2%s.pdf' %(date1+date2))

#t=np.arange(-12,-9,0.01)
#plt.figure()
#plt.xlabel('bias current I $[\mu A/cm^2]$')
#plt.ylabel('potential barrier')
#plt.plot(t,barrier(t,params2[0],params2[1]),'y')
#plt.plot(t,barrier(t,params2[2],params2[3]),'y')
#plt.plot(t,2*barrier(t,params2[0],params2[1]),'y')
#plt.plot(t,2*barrier(t,params2[2],params2[3]),'y')
#plt.plot(t,qbarrier(t,paramsq[0],paramsq[1],paramsq[2]),colorv[0])
#plt.plot(t,qbarrier(t,paramsq[3],paramsq[4],paramsq[5]),colorv[1])
#plt.plot(t,2*qbarrier(t,paramsq[0],paramsq[1],paramsq[2]),colorv[2])
#plt.plot(t,2*qbarrier(t,paramsq[3],paramsq[4],paramsq[5]),colorv[3])
#plt.plot(xs,params[1,:],colorv[0]+'o',label='burst to eq')
#plt.plot(xs,params[3,:],colorv[1]+'o',label='eq to burst')
#plt.plot(xs,2*params[1,:],colorv[2]+'o',label='2x burst to eq')
#plt.plot(xs,2*params[3,:],colorv[3]+'o',label='2x eq to burst')
#plt.legend()
#plt.savefig('barriercompdfnew%s.pdf' %(date1+date2))	
