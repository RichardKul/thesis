#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

#matplotlib.rcParams.update({'font.size': 22})

timefac=1000 #convert ms to s

date=['realfast16alcoarsewstf','realfast9acoarsetf','realfast23mtf','realfast19mtf']
D=[25,30,35,45]
l=len(D)


istart=1
ivalues=20

vec=np.zeros((l,20))
vecx=np.zeros((l,ivalues))
ii=0
offset=np.zeros(l,dtype=int)

for m in range(0,l):
	x=D[m]
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date[m],x,y),"r")
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
colorv=['g','y','b','r']
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ $[s^{-1}]$')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl],label='D=%.2f' %(D[n]/100))
plt.plot([0.16, 0.16], [5*10**(-1), 50000], color='black', linestyle='-',label='$I_{crit}$')
plt.plot([-0.01, -0.01], [5*10**(-1), 50000], color='black', linestyle='-')
plt.legend()
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [1,0,2,3]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('dneur25critsp%s.pdf' %(date[0]+date[1]))

vec=np.zeros((l,20))
vecx=np.zeros((l,ivalues))
ii=0
offset=np.zeros(l,dtype=int)

for m in range(0,l):
	x=D[m]
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date[m],x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,20-offset[ii]):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

plt.figure()
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor F')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl],label='D=%.2f' %(D[n]/100))
plt.plot([0.16, 0.16], [10**(-2), 30000], color='black', linestyle='-',label='$I_{crit}$')
plt.plot([-0.01, -0.01], [10**(-2), 30000], color='black', linestyle='-')
plt.legend()
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [1,0,2,3]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('fneur25critsp%s.pdf' %(date[0]+date[1]))

vec=np.zeros((l,20))
vecx=np.zeros((l,ivalues))
ii=0
offset=np.zeros(l,dtype=int)

for m in range(0,l):
	x=D[m]
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date[m],x,y),"r")
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

N=50000000
dt=0.0005
T=N*dt
file=open('/home/richard/mastergit/NetBeansProjects/detmodel/countI9a.txt',"r")
col=[]
for k in file:
	row=k.split()
	col.append(float(row[1]))
cola=np.array(col)*timefac
xburst=np.arange(-0.20,0.31,0.01)

plt.figure()
#colorv=['y','g','b','r','c']
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('average firing rate <v>$[s^{-1}]$')
#plt.yscale('log')
#plt.xscale('log')
plt.plot(xburst[12:51],cola[12:51]/T,label='spiking firing rate $v_0$',color='black')
for n in range(0,l):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl],label='D=%.2f' %(D[n]/100))
#plt.plot([0.165, 0.165], [5*10**(-1), 50000], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot([-0.022, -0.022], [5*10**(-1), 50000], color='black', linestyle='-')
#plt.xlim(-0.08,0.3)
plt.legend()
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [1,0,2,3]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('gneur25critsp%s.pdf' %(date[0]+date[1]))

