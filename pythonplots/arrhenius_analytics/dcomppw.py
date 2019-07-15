#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def deff(rp,rm,v):
	return (v**2*rp*rm)/((rp+rm)**3)

first=5
total=11
deffs=4

btoeq=np.zeros((deffs,total))
eqtob=np.zeros((deffs,total))

for k2 in range(first,first+total):
	x=[]
	ratefile = open('rate16m%d.txt' %k2,'r')
	for k4 in ratefile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	for k in range(0,deffs):
		btoeq[k][k2-first]=1/ax[k]
		eqtob[k][k2-first]=1/ax[k+deffs]

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

plt.plot(xfit,deff(btoeq[2-2,:],eqtob[2-2,:],v),'y')#,label='D=2,theory' )
plt.plot(xfit,deff(btoeq[3-2,:],eqtob[3-2,:],v),'g')#,label='D=3,theory' )
plt.plot(xfit,deff(btoeq[4-2,:],eqtob[4-2,:],v),'b')#,label='D=4,theory' )
plt.plot(xfit,deff(btoeq[5-2,:],eqtob[5-2,:],v),'r')#,label='D=5,theory' )
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xso,veco[0,:],'yo',label='D=2')
plt.plot(xs,vec[0,:],'go',label='D=3')
plt.plot(xs,vec[1,:],'bo',label='D=4')
plt.plot(xs,vec[2,:],'ro',label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('dcomppw16m.pdf')
