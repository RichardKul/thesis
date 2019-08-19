#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
Dvar=[30]
lvar=len(Dvar)
yvar=[4,13,3]
yvalues=len(yvar)

date1='realfast19jjem2st'
date0='realfast11jjem2sh'
date='realfast11jjem2st'
datefull=date0+date+date1
D1=[45]
l1=len(D1)
D0=[35]
l0=len(D0)
D=[40,50]
l=len(D)
ltot=l0+l+lvar+l1
vec=np.zeros((ltot,20))
ii=0
for x in Dvar:
	colvar,colxvar=[],[]
	for kk in range(0,yvalues):
		ystart=1
		jj=0
		while jj < kk:
			ystart=ystart+yvar[jj]
			jj=jj+1
		for y in range(ystart,ystart+yvar[kk]):
			file=open('/home/richard/outhome/d%s%d%d.txt' % (datevar[kk],x,y),"r")
			for k in file:
				row=k.split()
				colxvar.append(float(row[0]))
				colvar.append(float(row[1]))
	colxavar=np.array(colxvar)
	colavar=np.array(colvar)
	for z in range(0,20):
		vec[ii][z]=colavar[z]
	ii=ii+1

for x in D1:
	col1,colx1=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,20):
		vec[ii][z]=cola1[z]
	ii=ii+1

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
for n in range(0,lvar):
	plt.plot(colxavar,vec[n,:],label='D=%s' %Dvar[n])
for n in range(lvar,l1+lvar):
	plt.plot(colxa1,vec[n,:],label='D=%s' %D1[n-lvar])
for n in range(l1+lvar,l1+lvar+l0):
	plt.plot(colxa0,vec[n,:],label='D=%s' %D0[n-lvar-l1])
for n in range(l1+lvar+l0,ltot):
	plt.plot(colxa,vec[n,:],label='D=%s' %D[n-lvar-l1-l0])
handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,3,1,4]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.legend()
plt.savefig('dneur%s.pdf' %datefull)

vec=np.zeros((ltot,20))
ii=0

for x in Dvar:
	colvar,colxvar=[],[]
	for kk in range(0,yvalues):
		ystart=1
		jj=0
		while jj < kk:
			ystart=ystart+yvar[jj]
			jj=jj+1
		for y in range(ystart,ystart+yvar[kk]):
			file=open('/home/richard/outhome/f%s%d%d.txt' % (datevar[kk],x,y),"r")
			for k in file:
				row=k.split()
				colxvar.append(float(row[0]))
				colvar.append(float(row[1]))
	colxavar=np.array(colxvar)
	colavar=np.array(colvar)
	for z in range(0,20):
		vec[ii][z]=colavar[z]
	ii=ii+1

for x in D1:
	col1,colx1=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,20):
		vec[ii][z]=cola1[z]
	ii=ii+1

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
for n in range(0,lvar):
	plt.plot(colxavar,vec[n,:],label='D=%s' %Dvar[n])
for n in range(lvar,l1+lvar):
	plt.plot(colxa1,vec[n,:],label='D=%s' %D1[n-lvar])
for n in range(l1+lvar,l1+lvar+l0):
	plt.plot(colxa0,vec[n,:],label='D=%s' %D0[n-lvar-l1])
for n in range(l1+lvar+l0,ltot):
	plt.plot(colxa,vec[n,:],label='D=%s' %D[n-lvar-l1-l0])
#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')

handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,3,1,4]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('fneur%s.pdf' %datefull)

vec=np.zeros((ltot,20))
ii=0

for x in Dvar:
	colvar,colxvar=[],[]
	for kk in range(0,yvalues):
		ystart=1
		jj=0
		while jj < kk:
			ystart=ystart+yvar[jj]
			jj=jj+1
		for y in range(ystart,ystart+yvar[kk]):
			file=open('/home/richard/outhome/g%s%d%d.txt' % (datevar[kk],x,y),"r")
			for k in file:
				row=k.split()
				colxvar.append(float(row[0]))
				colvar.append(float(row[1]))
	colxavar=np.array(colxvar)
	colavar=np.array(colvar)
	for z in range(0,20):
		vec[ii][z]=colavar[z]
	ii=ii+1

for x in D1:
	col1,colx1=[],[]	
	for y in range(1,21):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,20):
		vec[ii][z]=cola1[z]
	ii=ii+1

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
for n in range(0,lvar):
	plt.plot(colxavar,vec[n,:],label='D=%s' %Dvar[n])
for n in range(lvar,l1+lvar):
	plt.plot(colxa1,vec[n,:],label='D=%s' %D1[n-lvar])
for n in range(l1+lvar,l1+lvar+l0):
	plt.plot(colxa0,vec[n,:],label='D=%s' %D0[n-lvar-l1])
for n in range(l1+lvar+l0,ltot):
	plt.plot(colxa,vec[n,:],label='D=%s' %D[n-lvar-l1-l0])
#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')

handles, labels = plt.gca().get_legend_handles_labels()
order = [0,2,3,1,4]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('gneur%s.pdf' %datefull)


