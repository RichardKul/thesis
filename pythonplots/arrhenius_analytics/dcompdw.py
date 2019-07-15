#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


def barrier(a,b,x):
	return a * x + b
def r(r0,U,D):
	return r0*np.exp(-U/D)
def deff(rp,up,rm,um,v,D):
	return (v**2*r(rp,up,D)*r(rm,um,D))/((r(rp,up,D)+r(rm,um,D))**3)

first=5
total=11
deffs=4

rp=np.zeros(total)
up=np.zeros(total)
rm=np.zeros(total)
um=np.zeros(total)

for k2 in range(first,first+total):
	x=[]
	paramfile = open('param16m%d.txt' %k2,'r')
	for k4 in paramfile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	rp[k2-first]=ax[0]
	up[k2-first]=ax[1]
	rm[k2-first]=ax[2]
	um[k2-first]=ax[3]
	


v=1.33




veco=np.zeros((1,16))
vec=np.zeros((3,20))
ii=0



for x in [20]:
	col=[]	
	for y in range(5,21):
		file=open('/home/richard/outhome/d16m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
	cola=np.array(col)
	for z in range(0,16):
		veco[0][z]=cola[z]

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
plt.xlim(0.5,3)
plt.ylim(5*10**(-2),2*10**3)
plt.yscale('log')

plt.plot(xfit,deff(rp,up,rm,um,v,2),'y')#,label='D=2,theory' )
plt.plot(xfit,deff(rp,up,rm,um,v,3),'g')#,label='D=3,theory' )
plt.plot(xfit,deff(rp,up,rm,um,v,4),'b')#,label='D=4,theory' )
plt.plot(xfit,deff(rp,up,rm,um,v,5),'r')#,label='D=5,theory' )
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xso,veco[0,:],'yo',label='D=2')
plt.plot(xs,vec[0,:],'go',label='D=3')
plt.plot(xs,vec[1,:],'bo',label='D=4')
plt.plot(xs,vec[2,:],'ro',label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('dcompdw16m.pdf')
