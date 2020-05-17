#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def comp(x,b,c,d,e):
	return b*x**3+c*x**2+d*x+e

timefac=1000
N=5000000
dt=0.005
T=N*dt
dvalues=6
ivalues2=51
count=np.zeros((dvalues,ivalues2))
ii=0
file=open('/home/richard/NetBeansProjects/detrinzel/countrinzelnoise5.txt',"r")
colx=[]
for k in file:
	col=[]
	row=k.split()
	colx.append(float(row[0]))
	for l in range(1,dvalues+1):	
		col.append(float(row[l]))
	cola=np.array(col)
	for l in range(0,dvalues):
		count[l][ii]=cola[l]
	ii=ii+1
colxa=np.array(colx)

date=['realrinzelrangelong26d1','realrinzelrange26d1','realrinzelrange26d1','realrinzelrangeshort26d1','realrinzelrangeshort26d1']
D=[200,250,300,400,500]
l=len(D)

istart=1
ivalues=10

vec=np.zeros((l,ivalues))
vecx=np.zeros((l,ivalues))
ii=0
offset=np.zeros(l,dtype=int)

for m in range(0,l):
	x=D[m]
	col1,colx1=[],[]	
	for y in range(istart,istart+ivalues):
		file=open('/home/richard/outhome/g%s%d%d.txt' % (date[m],x,y),"r")
		for k in file:
			row=k.split()
			colx1.append(float(row[0]))
			col1.append(float(row[1]))
		if len(col1)+offset[ii]<y-istart+1:
			offset[ii]=offset[ii]+1
	colxa1=np.array(colx1)
	cola1=np.array(col1)
	for z in range(0,ivalues-offset[ii]):
		vec[ii][z]=cola1[z]*timefac
		vecx[ii][z]=colxa1[z]
	ii=ii+1


colorv=['y','g','b','r','c','k']
t=np.arange(-21,-6,0.1)
plt.xlabel('bias current I')
plt.ylabel('firing rate')
plt.xlim(-16.2,-9)
#plt.ylim(1.2,1.4)
#plt.yscale('log')
for n in range(0,dvalues-1):
	nl=round(ivalues-offset[n])
	plt.plot(vecx[n,0:nl],vec[n,0:nl],label='D=%.0f' %(D[n]/10))
	plt.plot(colxa,count[n+1,:]/T*timefac,color='black')#,label='D=%s' %(D[n]/10))
	#plt.plot(t,comp(t,b[n],c[n],d[n],e[n]),colorv[n])
#plt.plot(t,comp(t,popt[0],popt[1]),label='linear appr %f %f'%(popt[0],popt[1]))
plt.legend()
plt.savefig('detmocountrinzelcompnew.pdf')


