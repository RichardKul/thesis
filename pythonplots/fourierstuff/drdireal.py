#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


vec=np.zeros((5,20))
gdiff=np.zeros((5,19))
deffmitt=np.zeros((5,19))
ii=0
istep=0.25

for x in [20]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g16m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(float(row[1]))
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	for n in range(0,19):
		gdiff[ii][n]=(cola[n+1]-cola[n])/istep
	ii=ii+1

for x in [25]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g29m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(float(row[1]))
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	for n in range(0,19):
		gdiff[ii][n]=(cola[n+1]-cola[n])/istep
	ii=ii+1

for x in [30,40,50]:
	col=[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g7m%d%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			col.append(float(row[1]))
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	for n in range(0,19):
		gdiff[ii][n]=(cola[n+1]-cola[n])/istep
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


xs=np.arange(-0.75,4.25,0.25)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.yscale('log')
plt.plot(xs,vec[0,:],label='D=2')
plt.plot(xs,vec[1,:],label='D=2.5')
plt.plot(xs,vec[2,:],label='D=3')
plt.plot(xs,vec[3,:],label='D=4')
plt.plot(xs,vec[4,:],label='D=5')
#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('gneurp29m5.pdf')

plt.figure()
xsmitt=np.arange(-0.625,4.125,0.25)
plt.ylabel('dr/dI')
#plt.yscale('log')
plt.plot(xsmitt,gdiff[0,:],label='D=2')
plt.plot(xsmitt,gdiff[1,:],label='D=2.5')
plt.plot(xsmitt,gdiff[2,:],label='D=3')
plt.plot(xsmitt,gdiff[3,:],label='D=4')
plt.plot(xsmitt,gdiff[4,:],label='D=5')

#plt.plot(xs,vec[4,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(xs,vec[6,:],label='D=4')

plt.legend()
plt.savefig('drdi29m5.pdf')
