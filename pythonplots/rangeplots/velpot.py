#!usr/bin/python3

import matplotlib

matplotlib.use("agg")

import numpy as np

import matplotlib.pyplot as plt

def velpot(v,F):
	return v**4/4-v**2/2-F*v
def lokex(arr,xmin,imin,imax,dx,p):
	minind=round((imin-xmin)/dx)
	maxind=round((imax-xmin)/dx)
	parta=arr[minind:maxind]
	if(p<=0):
		ind=np.where(parta==parta.min())
	else:
		ind=np.where(parta==parta.max())
	return xmin+(minind+ind[0][0])*dx

xmin=-1.8
xmax=1.8
dx=0.01
F=0.05

ymin=-0.4
ymax=0.2

t=np.arange(xmin,xmax,dx)
asy=velpot(t,F)
maxv=lokex(asy,xmin,-0.5,0.5,dx,1)
maxt=np.arange(maxv-0.5,maxv+0.5,dx)
minlv=lokex(asy,xmin,-1.5,-0.5,dx,-1)
minlt=np.arange(minlv,minlv+0.45,dx)
minrv=lokex(asy,xmin,0.5,1.5,dx,-1)
minrt=np.arange(minrv-0.6,minrv,dx)

plt.figure()

plt.ylabel('velocity potential U(v)',fontsize=12)
#plt.xlabel('velocity v')
plt.ylim(ymin,ymax)
plt.plot(t,velpot(t,F))
plt.plot(maxt,np.ones(100)*velpot(maxv,F),'black')
plt.plot(minlt,np.ones(45)*velpot(minlv,F),'black')
plt.plot(minrt,np.ones(60)*velpot(minrv,F),'black')
#plt.legend()
plt.arrow(maxv,velpot(maxv,F),0,ymin-velpot(maxv,F),head_length = 0.010, head_width = 0.025,color='black', length_includes_head = True,ls=':')
plt.arrow(minlv,velpot(minlv,F),0,ymin-velpot(minlv,F),head_length = 0.010, head_width = 0.025,color='black', length_includes_head = True,ls=':')
plt.arrow(minrv,velpot(minrv,F),0,ymin-velpot(minrv,F),head_length = 0.010, head_width = 0.025,color='black', length_includes_head = True,ls=':')
plt.arrow(maxv+0.48,velpot(maxv,F),0,velpot(minrv,F)-velpot(maxv,F),head_length = 0.010, head_width = 0.025,color='black', length_includes_head = True)
plt.arrow(maxv+0.48,velpot(minrv,F),0,velpot(maxv,F)-velpot(minrv,F),head_length = 0.010, head_width = 0.025,color='black', length_includes_head = True)
plt.arrow(maxv-0.5,velpot(maxv,F),0,velpot(minlv,F)-velpot(maxv,F),head_length = 0.010, head_width = 0.025,color='black', length_includes_head = True)
plt.arrow(maxv-0.5,velpot(minlv,F),0,velpot(maxv,F)-velpot(minlv,F),head_length = 0.010, head_width = 0.025,color='black', length_includes_head = True)
plt.text(maxv+0.55,(velpot(maxv,F)+velpot(minrv,F))/2, '$\Delta U_+$', fontdict=None,fontsize=15)
plt.text(maxv-1,(velpot(maxv,F)+velpot(minlv,F))/2, '$\Delta U_-$', fontdict=None,fontsize=15)
plt.text(maxv,ymin-0.03, '$v_0$', fontdict=None,fontsize=12)
plt.text(minrv,ymin-0.03, '$v_+$', fontdict=None,fontsize=12)
plt.text(minlv,ymin-0.03, '$v_-$', fontdict=None,fontsize=12)
plt.text(minrv+0.3,ymin-0.03, 'velocity $v$', fontdict=None,fontsize=12)
plt.xticks([])
plt.yticks([])
plt.savefig('velpot.pdf')

