#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

matplotlib.rcParams.update({'font.size': 16})

def qbarrier(x,a,b,c):
	return a*x**2+b*x+c
def func(x, a, b):
	return a * np.exp(-b * x)


date3='new'+'realfast11jjem2sh'
date2='new'+'realfast19jjem2st'

ivalues=20
l=5
D1=[]
D3=[]
D2=[35,45,40,50,30]
Dvar=[]
D=D1+D2+D3+Dvar
Da=np.array(D)
btoeq=np.zeros((l,ivalues))
eqtob=np.zeros((l,ivalues))
params=np.zeros((4,ivalues))
paramsq=np.zeros(6)
paramsqrate=np.zeros(6)
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

av=0.1/12
v0=0.47
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

istart=1	
xnew=np.arange(-5+istart,-5+istart+ivalues)*0.02

#barrier fit
popt,pcov = curve_fit(qbarrier, xnew, params[1,:])
paramsq[0]=popt[0]
paramsq[1]=popt[1]
paramsq[2]=popt[2]
popt,pcov = curve_fit(qbarrier, xnew, params[3,:])
paramsq[3]=popt[0]
paramsq[4]=popt[1]
paramsq[5]=popt[2]

#prefactor fit
popt,pcov = curve_fit(qbarrier, xnew, params[0,:])
paramsqrate[0]=popt[0]
paramsqrate[1]=popt[1]
paramsqrate[2]=popt[2]
popt,pcov = curve_fit(qbarrier, xnew, params[2,:])
paramsqrate[3]=popt[0]
paramsqrate[4]=popt[1]
paramsqrate[5]=popt[2]

t1=np.arange(-0.1,0.32,0.01)
plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('potential barrier')
plt.plot(xnew,params[1,:],'go',label='run to eq')
plt.plot(xnew,params[3,:],'ro',label='eq to run')
plt.plot(xnew,params[1,:],'g')
plt.plot(xnew,params[3,:],'r')
#plt.plot(t1,qbarrier(t1,paramsq[0],paramsq[1],paramsq[2]),'g')
#plt.plot(t1,qbarrier(t1,paramsq[3],paramsq[4],paramsq[5]),'r')
plt.legend()
plt.tight_layout()
plt.savefig('barriereal4linebig.pdf')

plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('rate prefactor')
plt.plot(xnew,params[0,:],'go',label='run to eq')
plt.plot(xnew,params[2,:],'ro',label='eq to run')
plt.plot(xnew,params[0,:],'g')
plt.plot(xnew,params[2,:],'r')
#plt.plot(t1,qbarrier(t1,paramsqrate[0],paramsqrate[1],paramsqrate[2]),'g')
#plt.plot(t1,qbarrier(t1,paramsqrate[3],paramsqrate[4],paramsqrate[5]),'r')
plt.legend()
plt.tight_layout()
plt.savefig('ratereal4linebig.pdf')

plt.figure()
plt.xlabel('bias current I $[\mu A/cm^2]$')
plt.ylabel('potential barrier')
plt.plot(xnew,params[1,:],'g^',label='run to eq')
plt.plot(xnew,params[3,:],'r^',label='eq to run')
plt.plot(xnew,2*params[1,:],'gs',label='2x run to eq')
plt.plot(xnew,2*params[3,:],'rs',label='2x eq to run')
plt.plot(xnew,params[1,:],'g')
plt.plot(xnew,params[3,:],'r')
plt.plot(xnew,2*params[1,:],'g')
plt.plot(xnew,2*params[3,:],'r')
plt.plot([0.165, 0.165], [0,9], color='black', linestyle='-',label='$I_{crit}$')
plt.plot([-0.022, -0.022], [0,9], color='black', linestyle='-')
plt.legend()
plt.tight_layout()
plt.savefig('barriercomprealfit4linecritbig.pdf')


