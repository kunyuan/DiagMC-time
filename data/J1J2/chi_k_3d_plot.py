#!/usr/bin/env python
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import read_data

def get_r(x,y,L):
    r = y*L+x
    return r

L = 16

Quans=["ChiK"]
key="ChiK"

Files = []
Files.append(read_data.read_array("1.65_mf/1.65_quantities.dat", Quans))

X = np.arange(0, L)
Y = np.arange(0, L)


Chi=[[0 for i in range(L)] for i in range(L)]
for x in X:
    for y in Y:
        Chi[x][y]=Files[0][key][0][get_r(x,y,L)].real

X, Y = np.meshgrid(X, Y)
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Chi, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False)

ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("Susceptibility")

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=L, aspect=5)

plt.savefig("3d_L16_1.50_Chi_k.pdf")
plt.show()
