#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def deff(rp,rm,v):
	return (v**2*rp*rm)/((rp+rm)**3)
def func(x, a, b):
	return a * np.exp(-b * x)
def deffana(r0p,r0m,ap,am,bp,bm,t,D,v):
	return (v**2*r(r0p,ap,bp,t,D)*r(r0m,am,bm,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**3)
def fanoana(r0p,r0m,ap,am,bp,bm,t,D,v):
	return (2*v*r(r0p,ap,bp,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**2)
def vfunc(x):
	return 1.204+174/(9.8*500)*x
first=5
total=11
deffs=3
epsilon = 0.00001

btoeq=np.zeros((deffs,total))
eqtob=np.zeros((deffs,total))
params=np.zeros((4,total))
paramsav=np.zeros((4,total))
params2=np.zeros(4)

for k2 in range(first,first+total):
	x=[]
	ratefile = open('rate29m%d.txt' %k2,'r')
	for k4 in ratefile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	for k in range(0,deffs):
		btoeq[k][k2-first]=1/ax[k]
		eqtob[k][k2-first]=1/ax[k+deffs]

v=1.35
xs=[1/3,1/4,1/5]
for k2 in range(first,first+total):
	popt,pcov = curve_fit(func, xs, btoeq[:,k2-first])
	params[0][k2-first]=popt[0]
	params[1][k2-first]=popt[1]
	popt,pcov = curve_fit(func, xs, eqtob[:,k2-first])
	params[2][k2-first]=popt[0]
	params[3][k2-first]=popt[1]
rbte=np.mean(params[0,:])
retb=np.mean(params[2,:])
for k2 in range(first,first+total):
	popt,pcov = curve_fit(func, xs, btoeq[:,k2-first],bounds=((rbte-epsilon,-np.inf), (rbte+epsilon,np.inf)))
	paramsav[0][k2-first]=popt[0]
	paramsav[1][k2-first]=popt[1]
	popt,pcov = curve_fit(func, xs, eqtob[:,k2-first],bounds=((retb-epsilon,-np.inf), (retb+epsilon,np.inf)))
	paramsav[2][k2-first]=popt[0]
	paramsav[3][k2-first]=popt[1]	
xnew=np.arange(0.5,3.25,0.25)
popt,pcov = curve_fit(barrier, xnew, paramsav[1,:])
params2[0]=popt[0]
params2[1]=popt[1]
popt,pcov = curve_fit(barrier, xnew, paramsav[3,:])
params2[2]=popt[0]
params2[3]=popt[1]
eqfile2 = open('paramsdf16m.txt','w')
for k4 in range(0,2): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.write('%.6f\n'%rbte) 
for k4 in range(2,4): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.write('%.6f\n'%retb) 
eqfile2.close() 
	


veco=np.zeros((1,16))
vec=np.zeros((4,20))
ii=0



for x in [25]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/d29m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	ii=ii+1
for x in [30,40,50]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/d7m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	ii=ii+1

#for x in [40]:
#	col=[]	
#	for y in range(5,21):
#		file=open('/home/richard/outhome/g17a%d%d.txt' % (x,y),"r")
#		for k in file:
#			row=k.split()
#			col.append(row[1])
#	cola=np.array(col)
#	for z in range(0,16):
#		vec[ii][z]=cola[z]
#	ii=ii+1

xfit=np.arange(0.5,3.25,0.25)
xso=np.arange(0.25,4.25,0.25)
xs=np.arange(-0.75,4.25,0.25)
xsh=np.arange(1,4.2,0.2)
plt.xlabel('bias current I')
plt.ylabel('$D_{eff}$')
plt.xlim(0,3.5)
plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')

ap=ax[0]
am=ax[3]
bp=ax[1]
bm=ax[4]
r0p=rbte
r0m=retb

t=np.arange(0,4,0.01)
plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,2.5,vfunc(2.5)),'y')
plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,3,vfunc(3)),'g')
plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,4,vfunc(4)),'b')
plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,5,vfunc(5)),'r')
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xs,vec[0,:],'yo',label='D=2.5')
plt.plot(xs,vec[1,:],'go',label='D=3')
plt.plot(xs,vec[2,:],'bo',label='D=4')
plt.plot(xs,vec[3,:],'ro',label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('dcompdf29mvarv.pdf')

fano=np.zeros((4,20))
ii=0
for x in [25]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f29m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,20):
		fano[ii][z]=cola[z]
	ii=ii+1
for x in [30,40,50]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f7m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,20):
		fano[ii][z]=cola[z]
	ii=ii+1

plt.figure()
plt.xlabel('bias current I')
plt.ylabel('Fano factor')
plt.yscale('log')
t=np.arange(-1,4,0.01)
plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,2.5,vfunc(2.5)),'y')
plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,3,vfunc(3)),'g')
plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,4,vfunc(4)),'b')
plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,5,vfunc(5)),'r')
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xs,fano[0,:],'yo',label='D=2.5')
plt.plot(xs,fano[1,:],'go',label='D=3')
plt.plot(xs,fano[2,:],'bo',label='D=4')
plt.plot(xs,fano[3,:],'ro',label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('fcompdf29mvarv.pdf')


