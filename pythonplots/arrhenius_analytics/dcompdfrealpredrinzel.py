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
def fit(x,a,b):
	return a/x+b
def comp(x,b,c,d,e):
	return b*x**3+c*x**2+d*x+e

matplotlib.rcParams.update({'font.size': 18})
timefac=1000

D2=[20,25,30]
D3=[]
D1=[15]
Dvar=[]
D=D1+D2+D3+Dvar
Da=np.array(D)
l1=len(D1)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l1+l2+l3+lvar
date1='realrinzelrange26d1'
date2='realrinzelrange26d1'


istart=1
ivalues=10

xnew=np.arange(172+istart,172+istart+ivalues)*0.25


file=open('/home/richard/mastergit/pythonplots/arrhenius_analytics/detmocountparamparam2.txt',"r")
col,colx=[],[]
for k in file:
	row=k.split()
	colx.append(float(row[0]))
	col.append(float(row[1]))
cola=np.array(col)
colxa=np.array(colx)



ivalues=10

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


#xs=np.arange(-5+istart,-5+istart+ivalues)*0.02
pv=np.arange(0,ivalues)

plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ [$s^{-1}$]')
#plt.xlim(0,3.5)
#plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')
#plt.xlim(-0.08,0.2)
Dvec=np.array([2,3,4,5])
#colorv=['g','y','b','r','c']
colorv=[ '#1f77b4', '#ff7f0e', '#2ca02c','#d62728','#9467bd']
#t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#for n in range(0,l):
#	plt.plot(t,deffana(axx[2],axx[5],axx[0],axx[3],axx[1],axx[4],t,Da[n]/100,av,v0),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot((pv-20.25)*0.8,deffred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],comp((pv-20.25)*0.8,fit(Dvec[n],colxa[0],cola[0]),fit(Dvec[n],colxa[1],cola[1]),fit(Dvec[n],colxa[2],cola[2]),fit(Dvec[n],colxa[3],cola[3])),Dvec[n])*timefac,colorv[n],label='D=%.0f'%(Dvec[n]))
	

plt.plot([-10.8, -10.8], [10**(-26), 10**22], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')


#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1,4]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='upper right', bbox_to_anchor=(0.65, 0.55))
plt.legend()#(loc='upper right', bbox_to_anchor=(0.8, 0.67))
plt.tight_layout()
plt.savefig('dcompdfpwnewpred4%s.pdf' %(date1+date2))


plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('firing rate r [$s^{-1}$]')
#plt.yscale('log')

#colorv=['b','r','g','y','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,vana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot((pv-20.25)*0.8,vred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],comp((pv-20.25)*0.8,fit(Dvec[n],colxa[0],cola[0]),fit(Dvec[n],colxa[1],cola[1]),fit(Dvec[n],colxa[2],cola[2]),fit(Dvec[n],colxa[3],cola[3])),Dvec[n])*timefac,colorv[n],label='D=%.0f'%(Dvec[n]))

#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right')
plt.legend()
plt.tight_layout()
plt.savefig('gcompdfpwnewpred4%s.pdf' %(date1+date2))


plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor F')
plt.yscale('log')
#plt.xlim(-0.08,0.2)
#colorv=['b','r','g','y','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot((pv-20.25)*0.8,fanored(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],comp((pv-20.25)*0.8,fit(Dvec[n],colxa[0],cola[0]),fit(Dvec[n],colxa[1],cola[1]),fit(Dvec[n],colxa[2],cola[2]),fit(Dvec[n],colxa[3],cola[3])),Dvec[n])*timefac,colorv[n],label='D=%.0f'%(Dvec[n]))
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.plot([-10.8, -10.8], [10**(-14), 10**25], color='black', linestyle='-',label='$I_{crit}$')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1,4]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower left')
plt.legend()
plt.tight_layout()
plt.savefig('fcompdfpwnewpred4%s.pdf' %(date1+date2))
