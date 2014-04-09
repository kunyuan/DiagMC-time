import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import read_data

#is2d = True
is2d = False
Beta = 0.9
N = 64
X = np.arange(0, Beta, Beta/N)
Y = np.arange(0, Beta, Beta/N)
X, Y = np.meshgrid(X, Y)
Gamt, dim, dim_name = read_data.read_array("./0_0.90_1_Gam_matrix.dat")
#print Gamt[0].diagonal().real.shape
#print X.shape
#print Y.shape

if is2d==True:
    tau = np.linspace(0, Beta, N, endpoint=False)
    fig = plt.figure()
    plt.plot(tau, Gamt[0].diagonal().real)
    plt.xlabel(dim_name[0])
    plt.ylabel("diag(Gamma)")
    plt.show()
    
else:
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(
        X, Y, Gamt[0].real, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    #csetz = ax.contour(X, Y, Gamt[0].real, zdir='Gamma', offset=0, cmap=cm.coolwarm)
    #csetx = ax.contour(X, Y, Gamt[0].real, zdir='x', offset=0,  cmap=cm.coolwarm)
    #csety = ax.contour(X, Y, Gamt[0].real, zdir='y', offset=Beta, cmap=cm.coolwarm)
    # ax.set_zlim(-1.01, 1.01)

    #surf = plt.plot(X, Gamt[0].real.diagonal())

    ax.set_xlabel(dim_name[0])
    ax.set_ylabel(dim_name[1])
    ax.set_zlabel("Gamma")

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=Beta, aspect=5)

    plt.show()
