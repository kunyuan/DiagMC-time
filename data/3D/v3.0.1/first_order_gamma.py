#!/usr/bin/env python
# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Is2D=True
#Is2D=False
N=64
Beta=1.0
t=np.arange(0,N,1)*Beta/N
G=(Beta-t)*np.exp(-1j*np.pi*t)
Gw=np.fft.fft(G)
Gw.resize(N,1)

Gamma=np.zeros((N,N),dtype=complex)
for ti in range(0,N):
	for tj in range(0,N):
		taui=ti*Beta/N
		tauj=tj*Beta/N
		Gamma[ti][tj]=(1.0+taui**2+tauj**2)*np.exp(-1j*np.pi*(taui+tauj))
Gammaw=np.fft.fft2(Gamma)
Gammaw=Gammaw+(N/Beta)**2
W=np.ones(N)
Ww=np.fft.fft(W)
GamOut=np.zeros((N,N),dtype=complex)
for wi in range(0,N):
    for wj in range(0,N):
        GamOut[wi][wj]=-Gw[wi][0]*Gw[wj][0]*Ww[0]*Gammaw[wi][wi]*Gammaw[wi][wj]*Gammaw[wj][wj]/N**7*Beta**6
Gamt=np.fft.ifft2(GamOut)
for ti in range(0,N):
    for tj in range(0,N):
        taui=ti*Beta/N
        tauj=tj*Beta/N
        Gamt[ti][tj]=Gamt[ti][tj]*np.exp(1j*np.pi*(taui+tauj))

if Is2D:
	fig=plt.figure()
	plt.plot(t,Gamt.diagonal().real,'r')
	plt.show()
else:
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	X = np.arange(0, Beta, 1.0/N)
	Y = np.arange(0, Beta, 1.0/N)
	X, Y = np.meshgrid(X, Y)

	surf = ax.plot_surface(X, Y, Gamt.real, rstride=1, cstride=1, cmap=cm.coolwarm,
	        linewidth=0, antialiased=False)
	# ax.set_zlim(-1.01, 1.01)

	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

	fig.colorbar(surf, shrink=0.5, aspect=5)

	plt.show()


