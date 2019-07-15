#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

veco=np.zeros((2,15))
vec=np.zeros((4,20))
ii=0

for x in [12]:
	col=[]	
	for y in range(5,6):
		file=open('/home/richard/outhome/d26a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(float(row[1]))
	for y in range(7,21):
		file=open('/home/richard/outhome/d26a%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(float(row[1]))
	cola=np.array(col)
	for z in range(0,15):
		veco[0][z]=cola[z]

for x in [25]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f29m%d%d.txt' % (x,y),"r")
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
		file=open('/home/richard/outhome/f7m%d%d.txt' % (x,y),"r")
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

xso=np.arange(0.5,4.25,0.25)
xs=np.arange(-0.75,4.25,0.25)
xsh=np.arange(1,4.2,0.2)
plt.xlabel('bias current I')
plt.ylabel('Fano factor')
plt.yscale('log')
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xs,vec[0,:],label='D=2.5')
plt.plot(xs,vec[1,:],label='D=3')
plt.plot(xs,vec[2,:],label='D=4')
plt.plot(xs,vec[3,:],label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('fneurp29mc.pdf')
