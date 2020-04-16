#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

#matplotlib.rcParams.update({'font.size': 22})

timefac=1000 #convert ms to s

date='realfast19mtf'
date0='realfast23mtf'
date1='realfast9acoarsetf'
D=[45]
D1=[30]
datefull=date+date1+date0
l1=len(D1)
D0=[25,35]
l0=len(D0)
l=len(D)
ltot=l0+l+l1
Dtot=D1+D0+D

istart=1
ivalues=20

vec=np.zeros((ltot,20))
vecx=np.zeros((ltot,ivalues))
ii=0
offset=np.zeros(ltot,dtype=int)

for x in D1:
	col1,colx1=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
		if len(col1)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola1[z]*timefac
		vecx[ii][z]=colxa1[z]
	ii=ii+1

for x in D0:
	col0,colx0=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date0,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
		if len(col0)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola0[z]*timefac
		vecx[ii][z]=colxa0[z]
	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola[z]*timefac
		vecx[ii][z]=colxa[z]
	ii=ii+1

plt.figure()
colorv=['y','g','b','r','c']
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ $[s^{-1}]$')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,ltot):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl],colorv[n],label='D=%.2f' %(Dtot[n]/100))
#plt.plot([0.165, 0.165], [5*10**(-1), 50000], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot([-0.022, -0.022], [5*10**(-1), 50000], color='black', linestyle='-')
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [1,0,2,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('dneur%s.pdf' %datefull)

vec=np.zeros((ltot,20))
vecx=np.zeros((ltot,ivalues))
ii=0
offset=np.zeros(ltot,dtype=int)

for x in D1:
	col1,colx1=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
		if len(col1)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola1[z]*timefac
		vecx[ii][z]=colxa1[z]
	ii=ii+1

for x in D0:
	col0,colx0=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date0,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
		if len(col0)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola0[z]*timefac
		vecx[ii][z]=colxa0[z]
	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola[z]*timefac
		vecx[ii][z]=colxa[z]
	ii=ii+1

plt.figure()
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor F')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,ltot):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl],colorv[n],label='D=%.2f' %(Dtot[n]/100))
#plt.plot([0.165, 0.165], [5*10**(-1), 50000], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot([-0.022, -0.022], [5*10**(-1), 50000], color='black', linestyle='-')
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [1,0,2,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('fneur%s.pdf' %datefull)

vec=np.zeros((ltot,20))
vecx=np.zeros((ltot,ivalues))
ii=0
offset=np.zeros(ltot,dtype=int)

for x in D1:
	col1,colx1=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
		if len(col1)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola1[z]*timefac
		vecx[ii][z]=colxa1[z]
	ii=ii+1

for x in D0:
	col0,colx0=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date0,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
		if len(col0)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola0[z]*timefac
		vecx[ii][z]=colxa0[z]
	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola[z]*timefac
		vecx[ii][z]=colxa[z]
	ii=ii+1

plt.figure()
colorv=['y','g','b','r','c']
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('average firing rate <v>$[s^{-1}]$')
#plt.yscale('log')
#plt.xscale('log')
for n in range(0,ltot):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl],colorv[n],label='D=%.2f' %(Dtot[n]/100))
#plt.plot([0.165, 0.165], [5*10**(-1), 50000], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot([-0.022, -0.022], [5*10**(-1), 50000], color='black', linestyle='-')
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [1,0,2,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('gneur%s.pdf' %datefull)

