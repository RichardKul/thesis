#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft


def barrier(a,b,x):
	return a * x + b
def r(r0,a,b,t,D):
	return r0*np.exp(-barrier(a,b,t)/D)
def drdi(t,r0,a,b):
	return (a*r0)/((1+np.exp(-barrier(a,b,t)))*(1+np.exp(barrier(a,b,t))))
def snr(r0p,r0m,ap,am,bp,bm,t,D):
	return (ap-am)**2/(D**2*(r(r0p,ap,bp,t,D)**(-1)+r(r0m,am,bm,t,D)**(-1)))
def snrtheo(t,r0,a,b,f):
	return drdi(t,r0,a,b)**2/(8*f**2)


#f=532812/1528
f=1

eqfile2 = open('/home/richard/mastergit/pythonplots/rangeplots/rateparams27a.txt','r')
rn,av,bv=[],[],[]
for k in eqfile2:
	row=k.split()
	rn.append(float(row[0]))
	av.append(float(row[1]))
	bv.append(float(row[2]))
rna=np.array(rn)
ava=np.array(av)
bva=np.array(bv)
rnaact=np.zeros(3)
avaact=np.zeros(3)
bvaact=np.zeros(3)
rnaact[0]=rna[4]
avaact[0]=ava[4]
bvaact[0]=bva[4]
rnaact[1]=rna[0]
avaact[1]=ava[0]
bvaact[1]=bva[0]
rnaact[2]=rna[1]
avaact[2]=ava[1]
bvaact[2]=bva[1]

bp = 1.74
ap = 5.64
r0p = 0.0075
bm = 3.15
am = -10.76
r0m = 0.012

date='realfast13aem2n4'
date1='realfast9aem2sh'
D=[30]
D1=[35,45]
Dtot=D+D1
l=len(D)+len(D1)

points=1000000
length=500000
ivalues=20
SNR=np.zeros((l,ivalues))
scale=np.zeros((l,ivalues))
ii=0

for c in D:
	for z in range(1,21):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date,c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date,c,z),"r")
		ll=0
		name,value=[],[]
		for k in param:
			row=k.split()
			lp=len(row)
			if ll<1:
				for jj in range(0,lp):
					name.append(row[jj])
			else:
				for kk in range(0,lp):
					value.append(float(row[kk]))
			ll=ll+1
		dt=value[name.index('dt')]
		N=value[name.index('N')]-value[name.index('Neq')]
		repetitions=value[name.index('repetitions')]
		epsilon=value[name.index('epsilon')]
		omega=value[name.index('omega')]
		T=N*repetitions*dt
		scale[ii][z-1]=epsilon**2*T
		omegaind=round(omega*T)		
		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1

for c1 in D1:
	for z in range(1,21):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date1,c1,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		ax=np.array(x)
		ay=np.array(y)
		param=open('/home/richard/outhome/param%s%d%d.txt' %(date1,c1,z),"r")
		ll=0
		name,value=[],[]
		for k in param:
			row=k.split()
			lp=len(row)
			if ll<1:
				for jj in range(0,lp):
					name.append(row[jj])
			else:
				for kk in range(0,lp):
					value.append(float(row[kk]))
			ll=ll+1
		dt=value[name.index('dt')]
		N=value[name.index('N')]-value[name.index('Neq')]
		repetitions=value[name.index('repetitions')]
		epsilon=value[name.index('epsilon')]
		omega=value[name.index('omega')]
		T=N*repetitions*dt
		scale[ii][z-1]=epsilon**2*T
		omegaind=round(omega*T)		
		SNR[ii][z-1]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
	ii=ii+1




D3=[35]
D2=[45]
Dvar=[30]
D=Dvar+D3+D2
Da=np.array(D)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l2+l3+lvar
date3='realfast11jjem2sh'
date2='realfast19jjem2st'
datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
yvar=[4,13,3]
yvalues=len(yvar)


ii=0

istart2=1
ivalues2=20

vec=np.zeros((l,ivalues2))
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

istart3=1
ivalues3=20

for x in D3:
	col3,colx3=[],[]	
	for y in range(istart3,istart3+ivalues3):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date3,x,y),"r")
		for k in file:
			row=k.split()
			colx3.append(float(row[0]))
			col3.append(float(row[1]))
	colxa3=np.array(colx3)
	cola3=np.array(col3)
	for z in range(0,ivalues3):
		vec[ii][z]=cola3[z]
	ii=ii+1

for x in D2:
	col2,colx2=[],[]	
	for y in range(istart2,istart2+ivalues2):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date2,x,y),"r")
		for k in file:
			row=k.split()
			colx2.append(float(row[0]))
			col2.append(float(row[1]))
	colxa2=np.array(colx2)
	cola2=np.array(col2)
	for z in range(0,ivalues2):
		vec[ii][z]=cola2[z]
	ii=ii+1





#files=open('/home/richard/NetBeansProjects/oup/xtraje6newwcav2.txt',"r")
#sx,sy=[],[]
#for ks in files:
#	rows=ks.split()
#	sx.append(float(rows[0]))
#	sy.append(float(rows[1]))
#sax=np.array(sx)
#say=np.array(sy)
#sax2=np.zeros(length) 
#say2=np.zeros(length)
#for l in range(0,length):
#	sax2[l]=sax[l]
#	say2[l]=say[l]
#plt.figure()
#plt.xlabel('time')
#plt.ylabel('position')
#plt.xlim(0,10)
#plt.plot(ax,ay)
#plt.plot(ax,aa)
#plt.savefig('oub.pdf')

#ys = fft(ay)
#for l in range(0,length):
#	S[l]=abs(ys[l])*abs(ys[l])/T

#omega=np.arange(0,length)*2*np.pi/T
plt.figure()
plt.xlabel('bias current')
plt.ylabel('SNR')

t=np.arange(-0.1,0.3,0.01)
xs=np.arange(-0.08,0.32,0.02)
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
colorv=['y','g','b','r','c']
for n in range(0,l):
	plt.plot(xs,(SNR[n,:]-1)/scale[n,:],colorv[n],label='D=%.2f' %(Dtot[n]*0.01))
for n in range(0,l):	
	plt.plot(xs,snrtheo(xs,rnaact[n],avaact[n],bvaact[n],f)/vec[n,:],colorv[n]+'o')
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
#handles, labels = plt.gca().get_legend_handles_labels()
#order = [0,2,1]
#plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
plt.legend()
plt.savefig('snrangerealanasp.pdf')
