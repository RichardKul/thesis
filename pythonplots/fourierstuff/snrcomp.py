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
def rprime(x,r0,a,D):
	return (1.204+174/(9.8*500)*x)*(-r0*a*np.exp(-a*x/D))/(D*(1+r0*np.exp(-a*x/D))**2)+(174/(9.8*500))/(1+r0*np.exp(-a*x/D))


first=5
total=11
deffs=3
epsilon = 0.00001

btoeq=np.zeros((deffs,total))
eqtob=np.zeros((deffs,total))
params=np.zeros((4,total))
paramsav=np.zeros((4,total))
params2=np.zeros(4)

file0=open('gpara2.txt',"r")
rvec,avec=[],[]
for k2 in file0:
	row=k2.split()
	rvec.append(float(row[0]))
	avec.append(float(row[1]))
rveca=np.array(rvec)
aveca=np.array(avec)

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
#plt.xlim(0,3.5)
#plt.ylim(5*10**(-2),2*10**3)
#plt.yscale('log')

#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xs,rprime(xs,rveca[0],aveca[0],2.5)**2/vec[0,:],'yo',label='D=2.5')
plt.plot(xs,rprime(xs,rveca[1],aveca[1],3)**2/vec[1,:],'go',label='D=3')
plt.plot(xs,rprime(xs,rveca[2],aveca[2],4)**2/vec[2,:],'bo',label='D=4')
plt.plot(xs,rprime(xs,rveca[3],aveca[3],5)**2/vec[3,:],'ro',label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('snrcompdf29m.pdf')




