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

D1=[30]
D2=[40,50]
D=D1+D2
Da=np.array(D)
l1=len(D1)
l2=len(D2)
l=l1+l2
date1='realfast6jjem2'
date2='realfast3jjem3'
istart=10
ivalues=6
epsilon = 0.00001

btoeq=np.zeros((l,ivalues))
eqtob=np.zeros((l,ivalues))
params=np.zeros((4,ivalues))
paramsav=np.zeros((4,ivalues))
params2=np.zeros(4)

for k2 in range(0,ivalues):
	x=[]
	ratefile = open('rate%s%d.txt' %(date1,k2),'r')
	for k4 in ratefile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	for k in range(0,l):
		btoeq[k][k2]=1/ax[k]
		eqtob[k][k2]=1/ax[k+l]

v=0.067
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
	
xnew=np.arange(-8+istart,-8+istart+ivalues)*0.025

popt,pcov = curve_fit(barrier, xnew, paramsav[1,:])
params2[0]=popt[0]
params2[1]=popt[1]
popt,pcov = curve_fit(barrier, xnew, paramsav[3,:])
params2[2]=popt[0]
params2[3]=popt[1]
eqfile2 = open('paramsdf%s.txt'%date1,'w')
for k4 in range(0,2): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.write('%.6f\n'%rbte) 
for k4 in range(2,4): 
	eqfile2.write('%.6f\n'%params2[k4]) 
eqfile2.write('%.6f\n'%retb) 
eqfile2.close() 
	




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

xs=np.arange(-8+istart1,-8+istart1+ivalues1)*0.025

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

colorv=['y','g','b','r']
t=np.arange(-0.2,0.3,0.01)
for n in range(0,l):
	plt.plot(t,deffana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot(xs,vec[n,:],colorv[n]+'o',label='D=%f' %(Da[n]/100))
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('dcompdf%s.pdf' %date1)

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


plt.figure()
plt.xlabel('bias current I')
plt.ylabel('Fano factor')
plt.yscale('log')

colorv=['y','g','b','r']
t=np.arange(-0.2,0.3,0.01)
for n in range(0,l):
	plt.plot(t,fanoana(rbte,retb,params2[1],params2[3],params2[0],params2[2],t,Da[n]/100,v),colorv[n])
#plt.plot(xsh,veco[0,:],label='D=1.2')
for n in range(0,l):
	plt.plot(xs,fano[n,:],colorv[n]+'o',label='D=%f' %(Da[n]/100))
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('fcompdf%s.pdf' %date1)


