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
def rate(r0p,r0m,ap,am,bp,bm,t,D,v):
	return v*r(r0m,am,bm,t,D)/(r(r0m,am,bm,t,D)+r(r0p,ap,bp,t,D))
def ratepar(r0p,r0m,ap,am,bp,bm,t):
	return 1.35*r(r0m,am,bm,t,4)/(r(r0m,am,bm,t,4)+r(r0p,ap,bp,t,4))
def sigmo(x,r0,a,D):
	return (1.204+174/(9.8*500)*x)/(1+r0*np.exp(-(a*x)/D))


first=5
total=11
deffs=3
epsilon = 0.00001

r0v=np.zeros(4)
av=np.zeros(4)
bv=np.zeros(4)
vv=np.zeros(4)
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
		file=open('/home/richard/outhome/g29m%d%d.txt' % (x,y),"r")
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
		file=open('/home/richard/outhome/g7m%d%d.txt' % (x,y),"r")
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
plt.ylabel('firing rate')
#plt.xlim(0,3.5)
#plt.ylim(5*10**(-2),2*10**3)
#plt.yscale('log')

ap=params2[1]
am=params2[3]
bp=params2[0]
bm=params2[2]
r0p=rbte
r0m=retb

t=np.arange(-1,4,0.01)
plt.plot(t,rate(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,2.5,v),'y')
plt.plot(t,rate(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,3,v),'g')
plt.plot(t,rate(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,4,v),'b')
plt.plot(t,rate(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,5,v),'r')
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xs,vec[0,:],'yo',label='D=2.5')
plt.plot(xs,vec[1,:],'go',label='D=3')
plt.plot(xs,vec[2,:],'bo',label='D=4')
plt.plot(xs,vec[3,:],'ro',label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('gcompdf29m.pdf')


#xnew=np.arange(-0.75,4.25,0.25)
#popt,pcov = curve_fit(ratepar, xnew, vec[2,:])
#r0p=popt[0]
#r0m=popt[1]
#ap=popt[2]
#am=popt[3]
#bp=popt[4]
#bm=popt[5]

#plt.figure()
#plt.xlabel('bias current I')
#plt.ylabel('firing rate')
#plt.plot(t,rate(r0p,r0m,ap,am,bp,bm,t,2.5,1.35),'y')
#plt.plot(t,rate(r0p,r0m,ap,am,bp,bm,t,3,1.35),'g')
#plt.plot(t,rate(r0p,r0m,ap,am,bp,bm,t,4,1.35),'b')
#plt.plot(t,rate(r0p,r0m,ap,am,bp,bm,t,5,1.35),'r')

#plt.plot(xs,vec[0,:],'yo',label='D=2.5')
#plt.plot(xs,vec[1,:],'go',label='D=3')
#plt.plot(xs,vec[2,:],'bo',label='D=4')
#plt.plot(xs,vec[3,:],'ro',label='D=5')

#plt.legend()
#plt.savefig('gfit29m.pdf')

xnew=np.arange(-0.75,4.25,0.25)
jj=0
for var in [2.5,3,4,5]:
	def sigmopar(x,r0,a):
		return (1.204+174/(9.8*500)*x)/(1+r0*np.exp(-(a*x)/var))
	popt,pcov = curve_fit(sigmopar, xnew, vec[jj,:])
	r0v[jj]=popt[0]
	av[jj]=popt[1]
	jj=jj+1
plt.figure()
plt.xlabel('bias current I')
plt.ylabel('firing rate')
plt.plot(t,sigmo(t,r0v[0],av[0],2.5),'y')
plt.plot(t,sigmo(t,r0v[1],av[1],3),'g')
plt.plot(t,sigmo(t,r0v[2],av[2],4),'b')
plt.plot(t,sigmo(t,r0v[3],av[3],5),'r')

plt.plot(xs,vec[0,:],'yo',label='D=2.5')
plt.plot(xs,vec[1,:],'go',label='D=3')
plt.plot(xs,vec[2,:],'bo',label='D=4')
plt.plot(xs,vec[3,:],'ro',label='D=5')
plt.legend()
plt.savefig('gfit29m.pdf')

eqfile3 = open('gpara2.txt','w')
for k3 in range(0,4): 
	eqfile3.write('%.6f %.6f %.6f %.6f \n'%(r0v[k3],av[k3],bv[k3],vv[k3])) 
eqfile3.close() 

