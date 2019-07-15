#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt



col,colx=[],[]	
for y in range(1,11):
	file=open('/home/richard/outhome/dgl%d.txt' % (y),"r")
	for k in file:
		row=k.split()
		col.append(float(row[1]))
		colx.append(float(row[0]))
cola=np.array(col)
colxa=np.array(colx)
plt.xlabel('gating variable $g_L$')
plt.ylabel('$D_{eff}$')
plt.yscale('log')
plt.plot(colxa,cola)

plt.savefig('dneurglvar.pdf')


