#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


jj=0
vecx=np.zeros((2,5,10000))
vecy=np.zeros((2,5,10000))
for x in ['ais','aism']:
	ii=0
	for y in range(1,6):
		colx,coly=[],[]
		file=open('/home/richard/outhome/deff%s%d.txt' % (x,y),"r")
		for k in file:
			row=k.split()
			colx.append(row[0])
			coly.append(row[1])
		colxa=np.array(colx)
		colya=np.array(coly)
		for z in range(0,10000):
			vecx[jj][ii][z]=colxa[z]
			vecy[jj][ii][z]=colya[z]
		ii=ii+1
	jj=jj+1

plt.xlabel('time')
plt.ylabel('$D_eff$')
plt.yscale('log')
plt.plot(vecx[0][0][:],vecy[0][0][:],'y',label='kT=0.033')
plt.plot(vecx[1][0][:],vecy[1][0][:],'y',label='kT=0.033, average')
plt.plot(vecx[0][1][:],vecy[0][1][:],'r',label='kT=0.043')
plt.plot(vecx[1][1][:],vecy[1][1][:],'r',label='kT=0.043, average')
plt.plot(vecx[0][2][:],vecy[0][2][:],'b',label='kT=0.056')
plt.plot(vecx[1][2][:],vecy[1][2][:],'b',label='kT=0.056, average')
plt.plot(vecx[0][3][:],vecy[0][3][:],'orange',label='kT=0.072')
plt.plot(vecx[1][3][:],vecy[1][3][:],'orange',label='kT=0.072, average')
plt.plot(vecx[0][4][:],vecy[0][4][:],'g',label='kT=0.094')
plt.plot(vecx[1][4][:],vecy[1][4][:],'g',label='kT=0.094, average')

#plt.legend()
plt.savefig('dmechav.pdf')


