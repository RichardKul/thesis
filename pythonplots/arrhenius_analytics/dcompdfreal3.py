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

matplotlib.rcParams.update({'font.size': 18})
timefac=1000

date=['realfast23mtf','realfast19mtf','realfast16alcoarsewstf','realfast9acoarsetf']
D=[35,45,25,30]
l=len(D)

ivalues=20

file=open('/home/richard/mastergit/NetBeansProjects/detmodel/countI9a.txt',"r")
col,colx=[],[]
for k in file:
    row=k.split()
    col.append(float(row[1]))
cola=np.array(col)
N=50000000
dt=0.0005
T=N*dt
vvec=np.zeros(ivalues)
for ll in range(0,ivalues):
	vvec[ll]=cola[2*ll+12]/T

#date10='realfast11jjem2sh'
#date20='realfast19jjem2st'
date10='realfast9acoarsetf'
date20='realfast23mtf'


rbtoeq=np.zeros(ivalues)
ubtoeq=np.zeros(ivalues)
reqtob=np.zeros(ivalues)
ueqtob=np.zeros(ivalues)

for k2 in range(0,ivalues):
    x=[]
    #ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/param%s%d.txt' %('new'+date10+'new'+date20,k2),'r')
    ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/param%s%d.txt' %(date10+date20,k2),'r')
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

ii=0

istart1=1
ivalues1=20

vec=np.zeros((l,ivalues1))

for m in range(0,l):
	x=D[m]
	col1,colx1=[],[]	
	for y in range(istart1,istart1+ivalues1):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date[m],x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues1):
		vec[ii][z]=cola1[z]*timefac
	ii=ii+1


xs=np.arange(-5+istart1,-5+istart1+ivalues1)*0.02
pv=np.arange(0,ivalues)

plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ [$s^{-1}$]')
#plt.xlim(0,3.5)
#plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')

colorv=[ '#2ca02c', '#d62728','#1f77b4', '#ff7f0e','#9467bd']
#colorv=['b','r','g','y','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#for n in range(0,l):
#	plt.plot(t,deffana(axx[2],axx[5],axx[0],axx[3],axx[1],axx[4],t,Da[n]/100,av,v0),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	if n!=2:
		plt.plot(xs,vec[n,:],'o',color=colorv[n],label='D=%.2f' %(D[n]/100))
	else:
		plt.plot(xs,vec[n,:],'o',color=colorv[n],label='D=%.2f$^*$' %(D[n]/100))
for n in range(0,l):
	plt.plot((pv-4)*0.02,deffred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],vvec[pv],D[n]/100)*timefac,colorv[n],label='D=%.2f'%(D[n]/100))
	
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')


handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='upper right', bbox_to_anchor=(0.65, 0.55))
plt.tight_layout()
plt.savefig('dcompdfpwnewbig2%s.pdf' %(date[0]+date[1]))

g=np.zeros((l,ivalues1))
ii=0

for m in range(0,l):
	x=D[m]
	col1,colx1=[],[]	
	for y in range(istart1,istart1+ivalues1):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date[m],x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues1):
		g[ii][z]=cola1[z]*timefac
	ii=ii+1

plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('firing rate r [$s^{-1}$]')
#plt.yscale('log')

colorv=[ '#2ca02c', '#d62728','#1f77b4', '#ff7f0e','#9467bd']
#colorv=['b','r','g','y','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,vana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	if n!=2:
		plt.plot(xs,g[n,:],'o',color=colorv[n],label='D=%.2f' %(D[n]/100))
	else:
		plt.plot(xs,g[n,:],'o',color=colorv[n],label='D=%.2f$^*$' %(D[n]/100))
for n in range(0,l):
	plt.plot((pv-4)*0.02,vred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],vvec[pv],D[n]/100)*timefac,colorv[n],label='D=%.2f'%(D[n]/100))

#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right')
plt.tight_layout()
plt.savefig('gcompdfpwnewbig2%s.pdf' %(date[0]+date[1]))


fano=np.zeros((l,ivalues1))
ii=0

for m in range(0,l):
	x=D[m]
	col1,colx1=[],[]	
	for y in range(istart1,istart1+ivalues1):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date[m],x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues1):
		fano[ii][z]=cola1[z]
	ii=ii+1

plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor F')
plt.yscale('log')

colorv=[ '#2ca02c', '#d62728','#1f77b4', '#ff7f0e','#9467bd']
#colorv=['b','r','g','y','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	if n!=2:
		plt.plot(xs,fano[n,:],'o',color=colorv[n],label='D=%.2f' %(D[n]/100))
	else:
		plt.plot(xs,fano[n,:],'o',color=colorv[n],label='D=%.2f$^*$' %(D[n]/100))
for n in range(0,l):
	plt.plot((pv-4)*0.02,fanored(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],vvec[pv],D[n]/100),colorv[n],label='D=%.2f'%(D[n]/100))
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')


handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower left')
plt.tight_layout()
plt.savefig('fcompdfpwnewbig2%s.pdf' %(date[0]+date[1]))

#t=np.arange(-0.1,0.3,0.001)
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
#plt.savefig('barriercompdfnewbig%s.pdf' %(date[0]+date[1]))
