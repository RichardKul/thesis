#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

fil1=open('/home/richard/outhome/g1s.txt',"r")
x1,y1=[],[]
for k in fil1:
    row=k.split()
    x1.append(float(row[0]))
    y1.append(float(row[1]))
x1s,y1s = zip(*sorted(zip(x1,y1)))

#fil2=open('/home/richard/outhome/g1.2mr.txt',"r")
#x2,y2=[],[]
#for k in fil2:
#    row=k.split()
#    x2.append(float(row[0]))
#    y2.append(float(row[1]))
#x2s,y2s = zip(*sorted(zip(x2,y2)))

#fil3=open('/home/richard/outhome/g1.5mr.txt',"r")
#x3,y3=[],[]
#for k in fil3:
#    row=k.split()
#    x3.append(float(row[0]))
#    y3.append(float(row[1]))
#x3s,y3s = zip(*sorted(zip(x3,y3)))

fil4=open('/home/richard/outhome/g2s.txt',"r")
x4,y4=[],[]
for k in fil4:
    row=k.split()
    x4.append(float(row[0]))
    y4.append(float(row[1]))
x4s,y4s = zip(*sorted(zip(x4,y4)))

fil5=open('/home/richard/outhome/g3s.txt',"r")
x5,y5=[],[]
for k in fil5:
    row=k.split()
    x5.append(float(row[0]))
    y5.append(float(row[1]))
x5s,y5s = zip(*sorted(zip(x5,y5)))

fil6=open('/home/richard/outhome/g4s.txt',"r")
x6,y6=[],[]
for k in fil6:
    row=k.split()
    x6.append(float(row[0]))
    y6.append(float(row[1]))
x6s,y6s = zip(*sorted(zip(x6,y6)))

#fil7=open('/home/richard/outhome/g5mr.txt',"r")
#x7,y7=[],[]
#for k in fil7:
#    row=k.split()
#    x7.append(float(row[0]))
#    y7.append(float(row[1]))
#x7s,y7s = zip(*sorted(zip(x7,y7)))

#fil8=open('/home/richard/outhome/g6mr.txt',"r")
#x8,y8=[],[]
#for k in fil8:
#    row=k.split()
#    x8.append(float(row[0]))
#    y8.append(float(row[1]))
#x8s,y8s = zip(*sorted(zip(x8,y8)))

plt.xlabel('bias current I')
plt.ylabel('firing rate')
plt.yscale('log')
plt.plot(x1s,y1s,label='D=1')
#plt.plot(x2s,y2s,label='D=1.2')
#plt.plot(x3s,y3s,label='D=1.5')
plt.plot(x4s,y4s,label='D=2')
plt.plot(x5s,y5s,label='D=3')
plt.plot(x6s,y6s,label='D=4')
#plt.plot(x7s,y7s,label='D=5')
#plt.plot(x8s,y8s,label='D=6')

plt.legend()
plt.savefig('gneurs.pdf')
