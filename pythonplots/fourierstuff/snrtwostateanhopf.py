#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft

from scipy.optimize import curve_fit

matplotlib.rcParams.update({'font.size': 14})

def func(x, a, b):
	return a * np.exp(-b * x)
def r(r0,U,D):
	return r0*np.exp(-U/D)
def snr(r0p,r0m,up,um,ups,ums,D,av,v0):
	return (r(r0m,um,D)*(r(r0p,up,D)+r(r0m,um,D))/(v0**2*r(r0p,up,D)))*(((ups-ums)*r(r0p,up,D)*v0)/(D*(r(r0p,up,D)+r(r0m,um,D)))+av)**2
def rprime(r0,r0s,Us,D):
	return (r0s/r0-Us/D)
def snrcor(r0p,r0m,up,um,ups,ums,D,av,v0,r0ps,r0ms):
	return ((av*r(r0m,um,D)*(r(r0m,um,D)+r(r0p,up,D))+v0*(rprime(r0m,r0ms,ums,D)-rprime(r0p,r0ps,ups,D))*r(r0p,up,D)*r(r0m,um,D))**2)/((r(r0p,up,D)+r(r0m,um,D))*v0**2*r(r0p,up,D)*r(r0m,um,D))

date3='realanhopf11flog'
date2='realanhopf19flog'

ivalues=14

D1=[20,25,30,35]
D2=[]
D3=[]
Dvar=[]
D=D1+D2+D3+Dvar
l=len(D1)
Da=np.array(D)
btoeq=np.zeros((l,ivalues))
eqtob=np.zeros((l,ivalues))
params=np.zeros((4,ivalues))
for k2 in range(0,ivalues):
	x=[]
	ratefile = open('/home/richard/mastergit/pythonplots/arrhenius_analytics/rate%s%d.txt' %(date3+date2,k2),'r')
	for k4 in ratefile:
		row=k4.split()
		x.append(float(row[0]))
	ax=np.array(x)
	for k in range(0,l):
		btoeq[k][k2]=1/ax[k]
		eqtob[k][k2]=1/ax[k+l]

xs=np.zeros(l)
for b in range(0,l):
	xs[b]=100/Da[b]
for k2 in range(0,ivalues):
	popt,pcov = curve_fit(func, xs, btoeq[:,k2])
	params[0][k2]=popt[0]
	params[1][k2]=popt[1]
	popt,pcov = curve_fit(func, xs, eqtob[:,k2])
	params[2][k2]=popt[0]
	params[3][k2]=popt[1]
rbte=np.mean(params[0,:])
retb=np.mean(params[2,:])

ups=np.zeros(ivalues-2)
ums=np.zeros(ivalues-2)
r0ms=np.zeros(ivalues-2)
r0ps=np.zeros(ivalues-2)

istart=4	
istep=0.25
xnew=np.arange(172+istart,172+istart+ivalues)*istep

for k3 in range(0,ivalues-2):
	ups[k3]=(params[1][k3+2]-params[1][k3])/(2*istep)
	ums[k3]=(params[3][k3+2]-params[3][k3])/(2*istep)
	r0ps[k3]=(params[0][k3+2]-params[0][k3])/(2*istep)
	r0ms[k3]=(params[2][k3+2]-params[2][k3])/(2*istep)

N0=50000000
dt0=0.005
T0=N0*dt0
file=open('/home/richard/mastergit/NetBeansProjects/detmodel/countanhopf.txt',"r")
col,colx=[],[]
for k in file:
	row=k.split()
	colx.append(float(row[0]))
	col.append(float(row[1]))
cola=np.array(col)/T0
colxa=np.array(colx)
lenx=len(colxa)
ratenew=np.zeros(ivalues)
dratenew=np.zeros(ivalues)
for ls in range(0,ivalues):
	x=xnew[ls]
	for m in range(0,lenx-1):		
		b=colxa[m+1]
		a=colxa[m]
		if x<=b and x>a:
			ratenew[ls]=cola[m]+((x-a)/(b-a))*(cola[m+1]-cola[m])
			dratenew[ls]=(cola[m+1]-cola[m])/(b-a)
			break

ll=0
lm=5
SNR=np.zeros((lm,ivalues))
snrfile=open('snrealanhopffile3.txt','r')
for s in snrfile:
	row=s.split()
	for t in range(istart-1-1,istart+ivalues-1-1):
		SNR[ll][t-istart+1+1]=float(row[t])	
	ll=ll+1

plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('Signal-to-noise ratio SNR')
Dtot=np.array([10,15,20,25,30])
l2=len(Dtot)
t=np.arange(-18,-8,0.1)
#xs=np.arange(-21.25+istart,-21.25+istart+ivalues)*0.8
#xs=np.arange(-20+istart,-20+istart+ivalues)*0.6
plt.yscale('log')
#plt.xscale('log')
#plt.xlim(4*10**(-3),5*10**3)
#plt.xlim(4*10**(-4),100)
colorv=['r','y','c','g','k','b'] # 6 colors
#colorv=['y','g','b'] # 3 colors
#colorv=['y','c','g','k','b'] # 5 colors
for n in range(1,l2):
	plt.plot(xnew,SNR[n,:],label='D=%.2f' %(Dtot[n]*0.01))
	#plt.plot(xnew[1:ivalues-1],snrcor(params[0,1:ivalues-1],params[2,1:ivalues-1],params[1,1:ivalues-1],params[3,1:ivalues-1],ups,ums,Dtot[n]*0.01,dratenew[1:ivalues-1],ratenew[1:ivalues-1],r0ps,r0ms)/8,colorv[n])#plt.plot(xnew,ratenew,colorv[n])
plt.plot([46.1, 46.1], [10**(-6), 10**(-2)], color='black', linestyle='-',label='$I_{crit}$')
plt.legend()
plt.tight_layout()
plt.savefig('snrtwostateanhopf7mnofitcrit.pdf')
#plt.savefig('snrinzelonly.pdf')
