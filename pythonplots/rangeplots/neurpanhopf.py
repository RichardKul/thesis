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

date='realanhopf26flog'
date1='realanhopf19flog'
date2='realrinzel15ninv0'
date3='realrinzel15ninv1'
D=[15]
D1=[20,25,30]
D2=[]
D3=[]
Dtot=D+D1+D2+D3
l=len(Dtot)
istart=4
ivalues=12
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
offset=np.zeros(l,dtype=int)
for x in D:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset[ii]):
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
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset[ii]):
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
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset[ii]):
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
		if len(col)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset[ii]):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

#xs=np.arange(-0.08,0.32,0.02)
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ $[s^{-1}]$')
plt.yscale('log')
#plt.xlim(44,47)
#plt.ylim(10**(-3),5*10**5)
#plt.xscale('log')
for n in range(0,l):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl]*timefac,label='D=%.2f' %(Dtot[n]/100))
#plt.plot([46.1, 46.1], [10**(-6), 10**3], color='black', linestyle='-')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

plt.legend()
plt.savefig('dneur3%s.pdf' %(date+date1))

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

offset2=np.zeros(l,dtype=int)
for x in D:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/f%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset2[ii]<y-istart+1:
			offset2[ii]=offset2[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset2[ii]):
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
		if len(col)+offset2[ii]<y-istart+1:
			offset2[ii]=offset2[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset2[ii]):
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
		if len(col)+offset2[ii]<y-istart+1:
			offset2[ii]=offset2[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset2[ii]):
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
		if len(col)+offset2[ii]<y-istart+1:
			offset2[ii]=offset2[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset2[ii]):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

plt.figure()

#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor')
plt.yscale('log')
#plt.xlim(44,47)
#plt.xscale('log')
for n in range(0,l):
	nl=round(ivalues-offset2[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl])#,label='D=%.2f' %(Dtot[n]/100))
#plt.plot([46.1, 46.1], [10**(-4), 10**5], color='black', linestyle='-')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

#plt.legend()
plt.savefig('fneur3%s.pdf' %(date+date1))

vec=np.zeros((l,ivalues))
vecx=np.zeros((l,ivalues))
ii=0

offset3=np.zeros(l,dtype=int)
for x in D:
	col,colx=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset3[ii]<y-istart+1:
			offset3[ii]=offset3[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset3[ii]):
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
		if len(col)+offset3[ii]<y-istart+1:
			offset3[ii]=offset3[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset3[ii]):
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
		if len(col)+offset3[ii]<y-istart+1:
			offset3[ii]=offset3[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset3[ii]):
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
		if len(col)+offset3[ii]<y-istart+1:
			offset3[ii]=offset3[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues-offset3[ii]):
		vec[ii][z]=cola[z]
		vecx[ii][z]=colxa[z]
	ii=ii+1

N=50000000
dt=0.005
T=N*dt
file=open('/home/richard/mastergit/NetBeansProjects/detmodel/countanhopf.txt',"r")
col,colx=[],[]
for k in file:
	row=k.split()
	colx.append(float(row[0]))
	col.append(float(row[1]))
colxa=np.array(colx)
cola=np.array(col)

plt.figure()
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I')
plt.ylabel('average firing rate <v> [$s^{-1}$]')
#plt.xlim(44,47)
#plt.yscale('log')
#plt.xscale('log')
plt.plot(colxa[10:38],cola[10:38]/T*timefac,label='running firing rate',color='black')
for n in range(0,l):	
	nl=round(ivalues-offset3[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl]*timefac)#,label='D=%.2f' %(Dtot[n]/100))
#plt.plot(colxa,vec[0,:],label='D=%.2f' %(D[0]/100))

#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')
plt.legend()
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [3,1,2,0]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])

plt.savefig('gneur3%s.pdf' %(date+date1))


