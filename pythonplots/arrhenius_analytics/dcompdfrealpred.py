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

matplotlib.rcParams.update({'font.size': 20})
timefac=1000
D1=[35]
D3=[25,30]
D2=[45]
Dvar=[]
D=D1+D2+D3+Dvar
Da=np.array(D)
l1=len(D1)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l1+l2+l3+lvar
date1='realfast11jjem2sh'
date3='realfast13aem2n4'
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
eqfile2 = open('paramsdfnewpred%s.txt'%(date1+date2),'w')
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

date10='realfast11jjem2sh'
date20='realfast19jjem2st'

ivalues=20

rbtoeq=np.zeros(ivalues)
ubtoeq=np.zeros(ivalues)
reqtob=np.zeros(ivalues)
ueqtob=np.zeros(ivalues)

for k2 in range(0,ivalues):
    x=[]
    ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/param%s%d.txt' %('new'+date10+'new'+date20,k2),'r')
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


xs=np.arange(-5+istart,-5+istart+ivalues)*0.02
pv=np.arange(0,ivalues)

plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ [$s^{-1}$]')
#plt.xlim(0,3.5)
#plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')
plt.xlim(-0.08,0.2)
Dvec=[2,3,4,5]
colorv=['g','y','b','r','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#for n in range(0,l):
#	plt.plot(t,deffana(axx[2],axx[5],axx[0],axx[3],axx[1],axx[4],t,Da[n]/100,av,v0),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot((pv-4)*0.02,deffred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],vvec[pv],Dvec[n]/100)*timefac,colorv[n],label='D=%.2f'%(Dvec[n]/100))
	
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.plot([0.165, 0.165], [10**(-35), 10**39], color='black', linestyle='-', label='$I_{crit}$')
plt.plot([-0.022, -0.022], [10**(-35), 10**39], color='black', linestyle='-')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1,4]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='upper right', bbox_to_anchor=(0.65, 0.55))
plt.legend(loc='upper right', bbox_to_anchor=(0.8, 0.67))
plt.tight_layout()
plt.savefig('dcompdfpwnewpred%s.pdf' %(date1+date2))


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
	plt.plot((pv-4)*0.02,vred(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],vvec[pv],Dvec[n]/100)*timefac,colorv[n],label='D=%.2f'%(Dvec[n]/100))

#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower right')
plt.legend()
plt.tight_layout()
plt.savefig('gcompdfpwnewpred%s.pdf' %(date1+date2))


plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor F')
plt.yscale('log')
plt.xlim(-0.08,0.2)
#colorv=['b','r','g','y','c']
t=np.arange(-0.1,0.3,0.01)
#for n in range(0,l):
#	plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot((pv-4)*0.02,fanored(rbtoeq[pv],ubtoeq[pv],reqtob[pv],ueqtob[pv],vvec[pv],Dvec[n]/100),colorv[n],label='D=%.2f'%(Dvec[n]/100))
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.plot([0.165, 0.165], [10**(-20), 10**40], color='black', linestyle='-', label='$I_{crit}$')
plt.plot([-0.022, -0.022], [10**(-20), 10**40], color='black', linestyle='-')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,3,0,1,4]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],loc='lower left')
plt.legend()
plt.tight_layout()
plt.savefig('fcompdfpwnewpred%s.pdf' %(date1+date2))
