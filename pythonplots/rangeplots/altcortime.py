#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


vec=np.zeros((3,20))
vecv=np.zeros((3,20))
ii=0


for x in [30]:
	col,colv=[],[]	
	for y in range(1,19):
		file=open('/home/richard/outhome/f1m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
		file=open('/home/richard/outhome/g1m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			colv.append(row[1])
	for y in range(18,19):
		file=open('/home/richard/outhome/f1m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
			col.append(row[1])
		file=open('/home/richard/outhome/g1m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			colv.append(row[1])
			colv.append(row[1])
	cola=np.array(col)
	colva=np.array(colv)
	for z in range(0,20):
		vec[ii][z]=cola[z]
		vecv[ii][z]=colva[z]
	ii=ii+1

for x in [40,50]:
	col,colv=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f1m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(row[1])
		file=open('/home/richard/outhome/g1m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			colv.append(row[1])
	cola=np.array(col)
	colva=np.array(colv)
	for z in range(0,20):
		vec[ii][z]=cola[z]
		vecv[ii][z]=colva[z]
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
plt.ylabel('Fano factor')
plt.yscale('log')
#plt.plot(xsh,veco[0,:],label='D=1.2')
plt.plot(xso,100*veco[1,:],label='D=2')
plt.plot(xs,vec[0,:],label='D=3')
plt.plot(xs,vec[1,:],label='D=4')
plt.plot(xs,vec[2,:],label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('fneurp1mcor.pdf')
