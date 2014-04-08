import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import read_data

Beta = 0.5
N = 64
fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(0, Beta, 0.50 / N)
Y = np.arange(0, Beta, 0.50 / N)
X, Y = np.meshgrid(X, Y)
Gamt, dim, dim_name = read_data.read_array("./0_0.50_1_Gam_matrix_new.dat")
print Gamt[0].real.shape
print X.shape
print Y.shape

surf = ax.plot_surface(
    X, Y, Gamt[0].real, rstride=1, cstride=1, cmap=cm.coolwarm,
    linewidth=0, antialiased=False)
# ax.set_zlim(-1.01, 1.01)

ax.set_xlabel(dim_name[0])
ax.set_ylabel(dim_name[1])
ax.set_zlabel("Gamma")

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
