#!/usr/bin/env python
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import read_data

Beta = 0.65
N = 64

#GamInt, dim_name = read_data.read_array("../../../project/0.50_quantities.dat")["Gamma"]
GamMC, dim_name = read_data.read_array("bare_L8_0.65_3_quantities.dat")["Gamma"]

X = np.arange(0, Beta, Beta/N)
Y = np.arange(0, Beta, Beta/N)
X, Y = np.meshgrid(X, Y)

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, GamMC.real, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False)

ax.set_xlabel(dim_name[0])
ax.set_ylabel(dim_name[1])
ax.set_zlabel("Gamma")

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=Beta, aspect=5)

plt.show()
