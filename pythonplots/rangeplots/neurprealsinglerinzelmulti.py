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

matplotlib.rcParams.update({'font.size': 14})
timefac=1000

date='realrinzel20ninv0'
date1='realrinzel15ninv0'
date2='realrinzel20ninv1'
date3='realrinzel15ninv1'
D=[]
D1=[200,300,500]
D2=[]
D3=[200,300,500]
Dtot=D+D1+D2+D3
l=len(D)+len(D1)+len(D2)+len(D3)
istart=2
ivalues=8
vecx=np.zeros((l,ivalues))
vec=np.zeros((l,ivalues))
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
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
for x in D1:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
for x in D2:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
for x in D3:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
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
for n in range(0,1):
	plt.plot(vecx[n,:],vec[n,:]*timefac,color='black',label='resting ICs')
for n in range(1,3):
	plt.plot(vecx[n,:],vec[n,:]*timefac,color='black')
for n in range(3,4):
	plt.plot(vecx[n,:],vec[n,:]*timefac,color='red',linestyle='--',label='spiking ICs')
for n in range(4,6):
	plt.plot(vecx[n,:],vec[n,:]*timefac,linestyle='--',color='red')
#,label='D=%.2f' %(Dtot[n]/10))
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [3,4,2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.tight_layout()
plt.savefig('dneursingle%s.pdf' %(date+date1))

vec=np.zeros((l,ivalues))
vecx=np.zeros((l,ivalues))
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
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D1:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D2:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
for x in D3:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
plt.figure()

#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor')
plt.yscale('log')
#plt.xscale('log')
#for n in range(0,l):
#	plt.plot(vecx[n,:],vec[n,:],label='D=%.2f' %(Dtot[n]/10))
for n in range(0,1):
	plt.plot(vecx[n,:],vec[n,:],color='black',label='resting ICs')
for n in range(1,3):
	plt.plot(vecx[n,:],vec[n,:],color='black')
for n in range(3,4):
	plt.plot(vecx[n,:],vec[n,:],color='red',label='spiking ICs')
for n in range(4,6):
	plt.plot(vecx[n,:],vec[n,:],color='red')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [3,4,2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('fneursingle%s.pdf' %(date+date1))

vec=np.zeros((l,ivalues))
vecx=np.zeros((l,ivalues))
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
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D1:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

for x in D2:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1
for x in D3:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues):
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
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('overall firing rate $<v> [s^{-1}]$')
#plt.yscale('log')
#plt.xscale('log')
#plt.plot(xburst,cola/T,label='measured bursting rate',color='black')
#for n in range(0,l):
#	plt.plot(vecx[n,:],vec[n,:],label='D=%.2f' %(Dtot[n]/10))
for n in range(0,1):
	plt.plot(vecx[n,:],vec[n,:]*timefac,color='black',label='resting ICs')
for n in range(1,3):
	plt.plot(vecx[n,:],vec[n,:]*timefac,color='black')
for n in range(3,4):
	plt.plot(vecx[n,:],vec[n,:]*timefac,linestyle='--',color='red',label='spiking ICs')
for n in range(4,6):
	plt.plot(vecx[n,:],vec[n,:]*timefac,linestyle='--',color='red')
#plt.plot(colxa,vec[0,:],label='D=%.2f' %(D[0]/100))

#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')
plt.legend()
plt.tight_layout()
plt.savefig('gneursingle%s.pdf' %(date+date1))


