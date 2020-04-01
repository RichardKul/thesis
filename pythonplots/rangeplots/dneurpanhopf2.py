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

date='realanhopf3mlogcorner2'
date1='realanhopf3mlogcorner'
date2='realanhopf21flogcorner'
date3='realanhopf3mlogcorner'
datevar=['realanhopf17mlogcorner2','realanhopf11mlogcorner2']

D=[]
D1=[10]
D2=[10]
D3=[]
Dvar=[10]
Dtot=D+D1+D2+D3+Dvar
l=len(Dtot)
istart=1
ivalues=10
vecx=np.zeros((l,ivalues))
vec=np.zeros((l,ivalues))
ii=0



offset=np.zeros(l,dtype=int)

yvar=[4,6]
yvalues=len(yvar)
ystart=1
for x in Dvar:
	colvar,colxvar=[],[]
	for kk in range(0,yvalues):
		for y in range(ystart,ystart+yvar[kk]):
			file=open('/home/richard/outhome/d%s%d%d.txt' % (datevar[kk],x,y),"r")
			for k in file:
				row=k.split()
				colxvar.append(float(row[0]))
				colvar.append(float(row[1]))
			if len(colvar)+offset[ii]<y+kk*yvar[kk]-istart+1:
				offset[ii]=offset[ii]+1
	colxavar=np.array(colxvar)
	colavar=np.array(colvar)
	for z in range(0,ivalues-offset[ii]):
		vec[ii][z]=colavar[z]
		vecx[ii][z]=colxavar[z]
	ii=ii+1
istart0=1
ivalues0=3
for x in D:
	col,colx=[],[]	
	for y in range(istart0,istart0+ivalues0):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx.append(float(row[0]))
			col.append(float(row[1]))
		if len(col)+offset[ii]<y-istart0+1:
			offset[ii]=offset[ii]+1
	colxa=np.array(colx)
	cola=np.array(col)
	for z in range(0,ivalues0-offset[ii]):
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
#plt.xscale('log')
for n in range(0,1):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl]*timefac,label='D=%.2f,kurz' %(Dtot[n]/100))#kurz
#for n in range(0,1):
#	plt.plot(vecx[n,0:3],vec[n,0:3]*timefac,label='D=%.2f' %(Dtot[n]/100))#kurz
#plt.xlim(46,47)
for n in range(1,2):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl]*timefac,label='D=%.2f,mittel' %(Dtot[n]/100))
for n in range(2,3):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl]*timefac,label='D=%.2f,lang' %(Dtot[n]/100))
#plt.plot([46.1, 46.1], [10**(-2), 10**2], color='black', linestyle='-',label='$I_{crit}$')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [2,0,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.legend()
plt.savefig('dneur%s.pdf' %(date+date1))


