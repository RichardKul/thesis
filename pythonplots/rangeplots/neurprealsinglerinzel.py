#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

#datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
#Dvar=np.array([30])
#lvar=len(Dvar)
#yvar=[4,13,3]
#yvalues=len(yvar)

timefac=1000

date='realrinzelrangelong26d1'
date1='realrinzelrange26d1'
date2='realrinzelrangeshort26d1'
D=[200]
D1=[250,300]
D2=[400,500]
Dtot=D+D1+D2
l=len(D)+len(D1)+len(D2)

vecx=np.zeros((l,10))
vec=np.zeros((l,10))
ii=0
#for x in Dvar:
#	colvar,colxvar=[],[]
#	for kk in range(0,yvalues):
#		ystart=1
#		jj=0
#		while jj < kk:
#			ystart=ystart+yvar[jj]
#			jj=jj+1
#		for y in range(ystart,ystart+yvar[kk]):
#			file=open('/home/richard/outhome/d%s%d%d.txt' % (datevar[kk],x,y),"r")
#			for k in file:
#				row=k.split()
#				colxvar.append(float(row[0]))
#				colvar.append(float(row[1]))
#	colxavar=np.array(colxvar)
#	colavar=np.array(colvar)
#	for z in range(0,20):
#		vec[ii][z]=colavar[z]
#	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
for x in D1:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
for x in D2:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

#xs=np.arange(-0.08,0.32,0.02)
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ $[s^{-1}]$')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l):
#	plt.plot(vecx[n,:],vec[n,:],label='D=%.2f,ic%i' %(Dtot[n]/10,round(n/3)))
	plt.plot(vecx[n,:],vec[n,:]*timefac,label='D=%.0f' %(Dtot[n]/10))
plt.plot([-10.8, -10.8], [10, 10**5], color='black', linestyle='-',label='$I_{crit}$')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [3,4,2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('dneursinglecrit%s.pdf' %(date+date1))

vec=np.zeros((l,10))
vecx=np.zeros((l,10))
ii=0

#for x in Dvar:
#	colvar,colxvar=[],[]
#	for kk in range(0,yvalues):
#		ystart=1
#		jj=0
#		while jj < kk:
#			ystart=ystart+yvar[jj]
#			jj=jj+1
#		for y in range(ystart,ystart+yvar[kk]):
#			file=open('/home/richard/outhome/f%s%d%d.txt' % (datevar[kk],x,y),"r")
#			for k in file:
#				row=k.split()
#				colxvar.append(float(row[0]))
#				colvar.append(float(row[1]))
#	colxavar=np.array(colxvar)
#	colavar=np.array(colvar)
#	for z in range(0,20):
#		vec[ii][z]=colavar[z]
#	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D1:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D2:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

plt.figure()

#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l):
	#plt.plot(vecx[n,:],vec[n,:],label='D=%.2f,ic%i' %(Dtot[n]/10,round(n/3)))
	plt.plot(vecx[n,:],vec[n,:],label='D=%.0f' %(Dtot[n]/10))
plt.plot([-10.8, -10.8], [0.1, 10**3], color='black', linestyle='-',label='$I_{crit}$')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [3,4,2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('fneursinglecrit%s.pdf' %(date+date1))

vec=np.zeros((l,10))
vecx=np.zeros((l,10))
ii=0

#for x in Dvar:
#	colvar,colxvar=[],[]
#	for kk in range(0,yvalues):
#		ystart=1
#		jj=0
#		while jj < kk:
#			ystart=ystart+yvar[jj]
#			jj=jj+1
#		for y in range(ystart,ystart+yvar[kk]):
#			file=open('/home/richard/outhome/g%s%d%d.txt' % (datevar[kk],x,y),"r")
#			for k in file:
#				row=k.split()
#				colxvar.append(float(row[0]))
#				colvar.append(float(row[1]))
#	colxavar=np.array(colxvar)
#	colavar=np.array(colvar)
#	for z in range(0,20):
#		vec[ii][z]=colavar[z]
#	ii=ii+1

for x in D:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D1:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D2:
	col,colx=[],[]	
	for y in range(1,11):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
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
cola=np.array(col)
xburst=np.arange(-0.20,0.31,0.01)

plt.figure()
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('firing rate [$s^{-1}$]')
#plt.yscale('log')
#plt.xscale('log')
#plt.plot(xburst,cola/T,label='measured bursting rate',color='black')
for n in range(0,l):
	#plt.plot(vecx[n,:],vec[n,:],label='D=%.2f,ic%i' %(Dtot[n]/10,round(n/3)))
	plt.plot(vecx[n,:],vec[n,:]*timefac,label='D=%.0f' %(Dtot[n]/10))
#plt.plot(colxa,vec[0,:],label='D=%.2f' %(D[0]/100))

#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')
plt.legend()
plt.savefig('gneursingle3%s.pdf' %(date+date1))


