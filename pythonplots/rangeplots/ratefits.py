#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def barrier(a,b,x):
	return a * x + b

def r(t,av,v0,a,b):
	return (av*t+v0)/(1+np.exp(-barrier(a,b,t)))

def ralt(t,v0,a,b):
	return r(t,0,v0,a,b)

N=50000000
dt=0.0005
T=N*dt
file=open('/home/richard/mastergit/NetBeansProjects/detmodel/countI9a.txt',"r")
col=[]
for k in file:
	row=k.split()
	col.append(float(row[1]))
cola=np.array(col)

D1=[35]
D3=[25,30]
D2=[45]
Dvar=[]
D=D1+D2+D3+Dvar
Da=np.array(D)
l1=len(D1)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l1+l2+l3+lvar
date1='realfast11jjem2sh'
date3='realfast13aem2n4'
date2='realfast19jjem2st'
datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
yvar=[4,13,3]
yvalues=len(yvar)

istart=1
ivalues=20
epsilon = 0.00001

g=np.zeros((l,ivalues))
xg=np.zeros((l,ivalues))
ii=0

istart1=1
ivalues1=20
for x in D1:
	col1,colx1=[],[]	
	for y in range(istart1,istart1+ivalues1):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date1,x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues1):
		g[ii][z]=cola1[z]
		xg[ii][z]=colxa1[z]
	ii=ii+1

istart2=1
ivalues2=20
for x in D2:
	col2,colx2=[],[]	
	for y in range(istart2,istart2+ivalues2):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx2.append(float(row[0]))
			col2.append(float(row[1]))
	colxa2=np.array(colx2)
	cola2=np.array(col2)
	for z in range(0,ivalues2):
		g[ii][z]=cola2[z]
		xg[ii][z]=colxa2[z]
	ii=ii+1

istart3=1
ivalues3=20
for x in D3:
	col3,colx3=[],[]	
	for y in range(istart3,istart3+ivalues3):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx3.append(float(row[0]))
			col3.append(float(row[1]))
	colxa3=np.array(colx3)
	cola3=np.array(col3)
	for z in range(0,ivalues3):
		g[ii][z]=cola3[z]
		xg[ii][z]=colxa3[z]
	ii=ii+1

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
		g[ii][z]=colavar[z]
		xg[ii][z]=colxavar[z]
	ii=ii+1

avn=np.zeros(l)
v0n=np.zeros(l)
av=np.zeros(l)
bv=np.zeros(l)

for kk in range(0,l):
	popt,pcov = curve_fit(r, xg[kk,:], g[kk,:],bounds=((0.01,0.06,-np.inf,-np.inf),(0.02,0.07,np.inf,np.inf)))
	avn[kk]=popt[0]
	v0n[kk]=popt[1]
	av[kk]=popt[2]
	bv[kk]=popt[3]

eqfile2 = open('rateparamsv11s.txt','w')
for k4 in range(0,l): 
	eqfile2.write('%.6f %.6f %.6f %.6f \n'%(avn[k4],v0n[k4],av[k4],bv[k4])) 
eqfile2.close() 

plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('average firing rate <v> $[10^3s^{-1}]$')
#plt.yscale('log')

colorv=['b','r','g','y','c']
t=np.arange(-0.1,0.31,0.01)
xburst=np.arange(-0.20,0.31,0.01)
#for n in range(0,l):
#	plt.plot(t,r(t,avn[n],v0n[n],av[n],bv[n]),colorv[n])
for n in range(0,l):
	plt.plot(xg[n,:],g[n,:],colorv[n],label='D=%.2f' %(Da[n]/100))
plt.plot(xburst,cola/T,label='running firing rate $v_0$',color='black')
plt.xlim(-0.1,0.3)
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,3,0,1,4]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.savefig('ganaburst23j2%s.pdf' %(date1+date2))
