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


date='realrinzel28o'

D=np.array([150,200])
l=len(D)

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
#			if x==300 and y==4:
#				colx.append(float(row[0])+0.8)
#				colx.append(float(row[0])+1.6)
#				col.append(float(row[1]))
#				col.append(float(row[1]))
#		colx.append(float(row[0]))
#		col.append(float(row[1]))
#	for y in range(5,7):
#		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
#		for k in file:
#			row=k.split()
#			colx.append(float(row[0]))
#			col.append(float(row[1]))
#		colx.append(float(row[0]))
#		col.append(float(row[1]))
#	for y in range(8,21):
#		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
#		for k in file:
#			row=k.split()
#			colx.append(float(row[0]))
#			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
	ii=ii+1

#xs=np.arange(-0.08,0.32,0.02)
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ $[10^3s^{-1}]$')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l):
	plt.plot(colxa,vec[n,:],label='D=%.2f' %(D[n]/10))
#plt.plot(colxa,vec[0,:],label='D=%.2f' %(D[0]/10))

plt.legend()
plt.savefig('dneursingle%s.pdf' %date)

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
#			if x==300 and y==4:
#				colx.append(float(row[0])+0.8)
#				colx.append(float(row[0])+1.6)
#				col.append(float(row[1]))
#				col.append(float(row[1]))
#		colx.append(float(row[0]))
#		col.append(float(row[1]))
#	for y in range(5,7):
#		file=open('/home/richard/outhome/f%s%d%d.txt' % (date,x,y),"r")
#		for k in file:
#			row=k.split()
#			colx.append(float(row[0]))
#			col.append(float(row[1]))
#		colx.append(float(row[0]))
#		col.append(float(row[1]))
#	for y in range(8,21):
#		file=open('/home/richard/outhome/f%s%d%d.txt' % (date,x,y),"r")
#		for k in file:
#			row=k.split()
#			colx.append(float(row[0]))
#			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
	ii=ii+1



plt.figure()

#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,l):
	plt.plot(colxa,vec[n,:],label='D=%.2f' %(D[n]/10))
#plt.plot(colxa,vec[0,:],label='D=%.2f' %(D[0]/10))
plt.legend()
plt.savefig('fneursingle%s.pdf' %date)

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
#			if x==300 and y==4:
#				colx.append(float(row[0])+0.8)
#				colx.append(float(row[0])+1.6)
#				col.append(float(row[1]))
#				col.append(float(row[1]))
#		colx.append(float(row[0]))
#		col.append(float(row[1]))
#	for y in range(5,7):
#		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
#		for k in file:
#			row=k.split()
#			colx.append(float(row[0]))
#			col.append(float(row[1]))
#		colx.append(float(row[0]))
#		col.append(float(row[1]))
#	for y in range(8,21):
#		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
#		for k in file:
#			row=k.split()
#			colx.append(float(row[0]))
#			col.append(float(row[1]))
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,10):
		vec[ii][z]=cola[z]
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
plt.ylabel('firing rate')
#plt.yscale('log')
#plt.xscale('log')
#plt.plot(xburst,cola/T,label='measured bursting rate',color='black')
for n in range(0,l):
	plt.plot(colxa,vec[n,:],label='D=%.2f' %(D[n]/10))
#plt.plot(colxa,vec[0,:],label='D=%.2f' %(D[0]/100))

#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')
plt.legend()
plt.savefig('gneursingle%s.pdf' %date)


