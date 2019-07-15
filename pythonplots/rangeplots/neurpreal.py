#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

date0='realfast6jjem2'
date='realfast3jjem3'
D0=[30]
l0=len(D0)
D=[40,50]
l=len(D)
ltot=l0+l
vec=np.zeros((ltot,20))
ii=0
for x in D0:
	col0,colx0=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date0,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,20):
		vec[ii][z]=cola0[z]
	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	ii=ii+1


#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('$D_{eff}$')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l0):
	plt.plot(colxa0,vec[n,:],label='D=%s' %D0[n])
for n in range(l0,ltot):
	plt.plot(colxa,vec[n,:],label='D=%s' %D[n-l0])

plt.legend()
plt.savefig('dneur%s.pdf' %date)

vec=np.zeros((ltot,20))
ii=0

for x in D0:
	col0,colx0=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date0,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,20):
		vec[ii][z]=cola0[z]
	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	ii=ii+1



plt.figure()
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('Fano factor')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l0):
	plt.plot(colxa0,vec[n,:],label='D=%s' %D0[n])
for n in range(l0,ltot):
	plt.plot(colxa,vec[n,:],label='D=%s' %D[n-l0])
#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')

plt.legend()
plt.savefig('fneur%s.pdf' %date)

vec=np.zeros((ltot,20))
ii=0

for x in D0:
	col0,colx0=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date0,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,20):
		vec[ii][z]=cola0[z]
	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,20):
		vec[ii][z]=cola[z]
	ii=ii+1



plt.figure()
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
#plt.yscale('log')
#plt.xscale('log')

for n in range(0,l0):
	plt.plot(colxa0,vec[n,:],label='D=%s' %D0[n])
for n in range(l0,ltot):
	plt.plot(colxa,vec[n,:],label='D=%s' %D[n-l0])
#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')

plt.legend()
plt.savefig('gneur%s.pdf' %date)


