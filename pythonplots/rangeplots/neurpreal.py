#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

#matplotlib.rcParams.update({'font.size': 22})

datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
Dvar=np.array([])
lvar=len(Dvar)
yvar=[4,13,3]
yvalues=len(yvar)

timefac=1000 #convert ms to s

date0='realfast11jjem2sh'
date='realfast19jjem2st'
date1='realfast13aem2n4'
D=np.array([45])
D1=np.array([30,25])
datefull=date+date1+date0
l1=len(D1)
D0=np.array([35])
l0=len(D0)
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
		vec[ii][z]=colavar[z]*timefac
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
		vec[ii][z]=cola1[z]*timefac
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
		vec[ii][z]=cola0[z]*timefac
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
		vec[ii][z]=cola[z]*timefac
	ii=ii+1

colorv=['y','g','b','r','c']
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('$D_{eff}$ $[s^{-1}]$')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,lvar):
	plt.plot(colxavar,vec[n,:],colorv[n],label='D=%.2f' %(Dvar[n]/100))
for n in range(lvar,l1+lvar):
	plt.plot(colxa1,vec[n,:],colorv[n],label='D=%.2f' %(D1[n-lvar]/100))
for n in range(l1+lvar,l1+lvar+l0):
	plt.plot(colxa0,vec[n,:],colorv[n],label='D=%.2f' %(D0[n-lvar-l1]/100))
for n in range(l1+lvar+l0,ltot):
	plt.plot(colxa,vec[n,:],colorv[n],label='D=%.2f' %(D[n-lvar-l1-l0]/100))
#plt.plot([0.165, 0.165], [5*10**(-1), 50000], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot([-0.022, -0.022], [5*10**(-1), 50000], color='black', linestyle='-')
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [1,0,2,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('dneur3sh%s.pdf' %datefull)

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
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Fano factor F')
plt.yscale('log')
#plt.xscale('log')
for n in range(0,lvar):
	plt.plot(colxavar,vec[n,:],colorv[n],label='D=%.2f' %(Dvar[n]/100))
for n in range(lvar,l1+lvar):
	plt.plot(colxa1,vec[n,:],colorv[n],label='D=%.2f' %(D1[n-lvar]/100))
for n in range(l1+lvar,l1+lvar+l0):
	plt.plot(colxa0,vec[n,:],colorv[n],label='D=%.2f' %(D0[n-lvar-l1]/100))
for n in range(l1+lvar+l0,ltot):
	plt.plot(colxa,vec[n,:],colorv[n],label='D=%.2f' %(D[n-lvar-l1-l0]/100))
#plt.plot([0.165, 0.165], [10**(-2), 10**4], color='black', linestyle='-', label='$I_{crit}$')
#plt.plot([-0.022, -0.022], [10**(-2), 10**4], color='black', linestyle='-')
#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [1,0,2,3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('fneur3sh%s.pdf' %datefull)

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
		vec[ii][z]=colavar[z]*timefac
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
		vec[ii][z]=cola1[z]*timefac
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
		vec[ii][z]=cola0[z]*timefac
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
		vec[ii][z]=cola[z]*timefac
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
#xs=np.arange(-0.75,4.25,0.25)
#xs=np.arange(0.25,2,0.25)
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('average firing rate <v> [$s^{-1}$]')
#plt.yscale('log')
#plt.xscale('log')
plt.plot(xburst,cola/T,label='running firing rate $v_0$',color='black')
for n in range(0,lvar):
	plt.plot(colxavar,vec[n,:],colorv[n],label='D=%.2f' %(Dvar[n]/100))
for n in range(lvar,l1+lvar):
	plt.plot(colxa1,vec[n,:],colorv[n],label='D=%.2f' %(D1[n-lvar]/100))
for n in range(l1+lvar,l1+lvar+l0):
	plt.plot(colxa0,vec[n,:],colorv[n],label='D=%.2f' %(D0[n-lvar-l1]/100))
for n in range(l1+lvar+l0,ltot):
	plt.plot(colxa,vec[n,:],colorv[n],label='D=%.2f' %(D[n-lvar-l1-l0]/100))
plt.xlim(-0.08,0.3)
#plt.plot(xs,vec[3,:],label='D=1.5')
#plt.plot(xs,vec[0,:],label='D=2')
#plt.plot(xs,vec[5,:],label='D=3')
#plt.plot(colxa,cola,label='D=3e-3')
#plt.legend()
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,1,3,4,0]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('gneursh3%s.pdf' %datefull)


