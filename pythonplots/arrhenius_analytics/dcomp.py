#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def deff(r0p,r0m,ap,am,bp,bm,t,D,v):
	return (v**2*r(r0p,ap,bp,t,D)*r(r0m,am,bm,t,D))/((r(r0p,ap,bp,t,D)+r(r0m,am,bm,t,D))**3)

file=open('/home/richard/NetBeansProjects/parambarrier16m.txt',"r")
x=[]
for k in file:
	row=k.split()
	x.append(float(row[0]))
ax=np.array(x)

t=np.arange(0,4,0.01)
ap=ax[0]
am=ax[3]
bp=ax[1]
bm=ax[4]
r0p=ax[2]
r0m=ax[5]
v=1.33




veco=np.zeros((2,16))
vec=np.zeros((3,20))
ii=0

for x in [12]:
	col=[]	
	for y in range(5,6):
		file=open('/home/richard/outhome/d26a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(float(row[1]))
			col.append(float(row[1]))
	for y in range(7,21):
		file=open('/home/richard/outhome/d26a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(float(row[1]))
	cola=np.array(col)
	for z in range(0,16):
		veco[0][z]=cola[z]

for x in [20]:
	col=[]	
	for y in range(5,21):
		file=open('/home/richard/outhome/d16m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,16):
		veco[1][z]=cola[z]

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

xso=np.arange(0.25,4.25,0.25)
xs=np.arange(-0.75,4.25,0.25)
xsh=np.arange(1,4.2,0.2)
plt.xlabel('bias current I')
plt.ylabel('$D_{eff}$')
plt.xlim(0.5,3)
plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')
plt.plot(t,deff(r0p,r0m,ap,am,bp,bm,t,2,v),'y')
plt.plot(t,deff(r0p,r0m,ap,am,bp,bm,t,3,v),'g')
plt.plot(t,deff(r0p,r0m,ap,am,bp,bm,t,4,v),'b')
plt.plot(t,deff(r0p,r0m,ap,am,bp,bm,t,5,v),'r')#,label='D=%d,theory' %k)
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xso,veco[1,:],'yo',label='D=2')
plt.plot(xs,vec[0,:],'go',label='D=3')
plt.plot(xs,vec[1,:],'bo',label='D=4')
plt.plot(xs,vec[2,:],'ro',label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('dcomp16mp.pdf')
