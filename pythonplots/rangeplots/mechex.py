#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt


vecx=np.zeros((5,30))
vecy=np.zeros((5,30))
ii=0
colx,coly=[],[]
for y in range(1,31):
	file=open('/home/richard/outhome/dex%dkt0.txt' % (y),"r")
	for k in file:
		row=k.split()
		colx.append(row[0])
		coly.append(row[1])
colxa=np.array(colx)
colya=np.array(coly)
for z in range(0,30):
	vecx[ii][z]=colxa[z]
	vecy[ii][z]=colya[z]
ii=ii+1
for x in range(1,5):
	colx,coly=[],[]
	for y in range(1,31):
		file=open('/home/richard/outhome/mechdex%dkt%d.txt' % (y,x),"r")
		for k in file:
			row=k.split()
			colx.append(row[0])
			coly.append(row[1])
	colxa=np.array(colx)
	colya=np.array(coly)
	for z in range(0,30):
		vecx[ii][z]=colxa[z]
		vecy[ii][z]=colya[z]
	ii=ii+1
	

plt.xlabel('bias force F')
plt.ylabel('$D_{eff}$')
plt.yscale('log')
plt.ylim(1,3*10**4)
plt.plot(vecx[0,:],vecy[0,:],label='kT=0.033')
plt.plot(vecx[1,:],vecy[1,:],label='kT=0.043')
plt.plot(vecx[2,:],vecy[2,:],label='kT=0.056')
plt.plot(vecx[3,:],vecy[3,:],label='kT=0.072')
plt.plot(vecx[4,:],vecy[4,:],label='kT=0.094')

plt.legend()
plt.savefig('dmechexaktnew.pdf')


