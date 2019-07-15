#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

fil1=open('/home/richard/outhome2/outhome/Gesch33sh2.txt',"r")
x1,y1=[],[]
for k in fil1:
    row=k.split()
    x1.append(float(row[0]))
    y1.append(float(row[1]))
x1s,y1s = zip(*sorted(zip(x1,y1)))

col=[]	
for y in range(1,31):
	file=open('/home/richard/outhome/mechg%d.txt' % (y),"r")
	for k in file:
		row=k.split()
		col.append(float(row[1]))
cola=np.array(col)
x1sa=np.arange(0.56,0.86,0.01)

fil2=open('/home/richard/outhome2/outhome/Gesch43sh2.txt',"r")
x2,y2=[],[]
for k in fil2:
    row=k.split()
    x2.append(float(row[0]))
    y2.append(float(row[1]))
x2s,y2s = zip(*sorted(zip(x2,y2)))

fil3=open('/home/richard/outhome2/outhome/Gesch56sh.txt',"r")
x3,y3=[],[]
for k in fil3:
    row=k.split()
    x3.append(float(row[0]))
    y3.append(float(row[1]))
x3s,y3s = zip(*sorted(zip(x3,y3)))

fil4=open('/home/richard/outhome2/outhome/Gesch72sh.txt',"r")
x4,y4=[],[]
for k in fil4:
    row=k.split()
    x4.append(float(row[0]))
    y4.append(float(row[1]))
x4s,y4s = zip(*sorted(zip(x4,y4)))

fil5=open('/home/richard/outhome2/outhome/Gesch94sh.txt',"r")
x5,y5=[],[]
for k in fil5:
    row=k.split()
    x5.append(float(row[0]))
    y5.append(float(row[1]))
x5s,y5s = zip(*sorted(zip(x5,y5)))

#xsa=np.array(xs)
#ysa=np.array(ys)
plt.xlabel('bias Force F')
plt.ylabel('velocity')
plt.plot(x1s,y1s,label='kT=0.033')
plt.plot(x2s,y2s,label='kT=0.043')
plt.plot(x3s,y3s,label='kT=0.056')
plt.plot(x4s,y4s,label='kT=0.072')
plt.plot(x5s,y5s,label='kT=0.094')
plt.legend()
plt.savefig('mechgsh.pdf')