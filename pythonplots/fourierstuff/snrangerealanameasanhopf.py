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
def snrmeas(drdr,f):
	return (drdr**2/(8*f**2)) 


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

#date='realfast13aem2n4'
#date1='realfast9aem2sh'
#D=[25,30]
#D1=[35,45]

date='realanhopf3fwsshort'
date1='realanhopf7flog'
date0='realrinzelrangeshort26d1'

D=[20,30]
D1=[10]
D0=[]
Dtot=D+D1+D0
l2=len(D)+len(D1)+len(D0)

#points=1000000
#length=500000
ivalues=20
istart=1
SNR=np.zeros((l2,ivalues))
scale=np.zeros((l2,ivalues))
xvec=np.zeros((l2,ivalues))
offset=np.zeros(l2,dtype=int)
ii=0

for c in D:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date,c,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		if len(x) == 0:
			offset[ii]=offset[ii]+1
			continue
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
		scale[ii][z-istart-offset[ii]]=epsilon**2*T
		omegaind=round(omega*T)	
		#if c==200 and z==9:
		#	omegaind=2105	
		#omegaind=141	
		SNR[ii][z-istart-offset[ii]]=ay[omegaind]/np.mean([ay[omegaind-1],ay[omegaind-2],ay[omegaind+1],ay[omegaind+2],ay[omegaind-3]])
		iv=open('/home/richard/outhome/d%s%d%d.txt' %(date,c,z),"r")
		for k in iv:
			row=k.split()
		xvec[ii][z-istart-offset[ii]]=float(row[0])
	ii=ii+1

for c1 in D1:
	for z in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/spike%s%d%d.txt' %(date1,c1,z),"r")
		x,y=[],[]
		for k in file:
			row=k.split()
			x.append(float(row[0]))
			y.append(float(row[1]))
		le=len(x)
		if le == 0:
			offset[ii]=offset[ii]+1
			continue
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
		scale[ii][z-istart-offset[ii]]=epsilon**2*T
		#omegaind=round(omega*T)
		omegaind=le-6		
		#omegaind=141
		SNR[ii][z-istart-offset[ii]]=2*ay[omegaind]/np.mean(ay[omegaind-3:omegaind]+ay[omegaind+1:omegaind+4])
		iv=open('/home/richard/outhome/d%s%d%d.txt' %(date1,c1,z),"r")
		for k in iv:
			row=k.split()
		xvec[ii][z-istart-offset[ii]]=float(row[0])
	ii=ii+1

snrealfile=open('snrealanhopffile.txt','w')
for n in range(0,l2):
	for m in range(0,ivalues):
		snrealfile.write('%.6f '%((SNR[n][m]-1)/scale[n][m]))
	snrealfile.write('\n')
snrealfile.close()

D0=[]
D3=[20,30,10]
D2=[]
Dvar=[]
D=D0+Dvar+D3+D2
Da=np.array(D)
l0=len(D0)
l2=len(D2)
l3=len(D3)
lvar=len(Dvar)
l=l2+l3+lvar+l0
date3='realanhopf22j'
date2='realfast19jjem2st'
datevar=['realfast11jjem2','realfast11jjem2sh','realfast11jjem2']
yvar=[4,13,3]
yvalues=len(yvar)


ii=0

istart0=1
ivalues0=20
istart2=1
ivalues2=20
vec=np.zeros((l,ivalues))

for x in D0:
	col0,colx0=[],[]	
	for y in range(istart0,istart0+ivalues0):
		file=open('/home/richard/outhome/d%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,ivalues0):
		vec[ii][z]=cola0[z]
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

istep=0.25

ii=0
gvec=np.zeros((l,ivalues))
drdr=np.zeros((l,ivalues-2))

for x in D0:
	col0,colx0=[],[]	
	for y in range(istart0,istart0+ivalues0):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date,x,y),"r")
		for k in file:
			row=k.split()
			colx0.append(float(row[0]))
			col0.append(float(row[1]))
	colxa0=np.array(colx0)
	cola0=np.array(col0)
	for z in range(0,ivalues0):
		gvec[ii][z]=cola0[z]
	for z in range(0,18):
		drdr[ii][z]=(cola0[z+2]-cola0[z])/(2*istep)
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
		gvec[ii][z]=colavar[z]
	for z in range(0,18):
		drdr[ii][z]=(colavar[z+2]-colavar[z])/(2*istep)
	ii=ii+1


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
		gvec[ii][z]=cola3[z]
	for z in range(0,18):
		drdr[ii][z]=(cola3[z+2]-cola3[z])/(2*istep)
	ii=ii+1

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
		gvec[ii][z]=cola2[z]
	for z in range(0,18):
		drdr[ii][z]=(cola2[z+2]-cola2[z])/(2*istep)
	ii=ii+1

plt.figure()
xsh=np.arange(43.5,48,0.25)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
plt.yscale('log')
for n in range(0,l):
	plt.plot(xsh,drdr[n,:],label='D=%s' %Dtot[n])
plt.legend()
plt.savefig('drdineurallanhopf.pdf')
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
xs=np.arange(43.25,48.25,0.25)

plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
colorv=['y','g','b','r']
for n in range(0,l-1):
	nl=round(ivalues-offset[n])
	plt.plot(xvec[n,0:nl],(SNR[n,0:nl]-1)/scale[n,0:nl],colorv[n]+'o',label='D=%.2f' %(Dtot[n]*0.1))
for n in range(l-1,l):
	nl=round(ivalues-offset[n])
	mult=np.ones(nl)
	for f2 in range(nl-7,nl):
		mult[f2]=10
	plt.plot(xvec[n,0:nl],mult*2*(SNR[n,0:nl]-1)/scale[n,0:nl],colorv[n]+'o',label='D=%.2f' %(Dtot[n]*0.1))
for n in range(0,l):	
	plt.plot(xsh,snrmeas(drdr[n,:],f)/vec[n,1:19],colorv[n])
#plt.plot([0.163, 0.163], [2*10**(-5), 6*10**(-2)], color='black', linestyle='-')
#plt.plot([-0.02, -0.02], [2*10**(-5), 6*10**(-2)], color='black', linestyle='-',label='$I_{crit}$')
#plt.plot(xs,SNR[2,:],label='D=3')
#plt.plot(xs,SNR[1,:],label='D=2.5')
handles, labels = plt.gca().get_legend_handles_labels()
order = [2,0,1]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
#plt.plot(sax2,say2/T2,label='e6')
#plt.legend()
plt.savefig('snrangerealanameasspallanhopf.pdf')
