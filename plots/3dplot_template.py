#!/usr/bin/env python
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import read_data

#is2d = True
is2d = False

Beta = 0.90
N = 64

#Quans = ["Gamma2"]
#Quans = ["Gamma"]
GamInt, dim_name = read_data.read_array("../project/0.90_Gam1.dat")["Gamma"]
#GamMC, dim_name = read_data.read_array("./../../data/conservation/bare_0.90_4_quantities.dat")["Gamma2"]
#GamMC, dim_name = read_data.read_array("./1_loop/0.90_quantities.dat")["Gamma"]
#GamMC = read_data.read_array("../0.90_quantities.dat", Quans)

if is2d is True:
    tau = np.arange(0, Beta, Beta/N)
    fig = plt.figure()
    #plt.plot(tau, GamMC.diagonal().real, 'r', 
    plt.plot(tau, GamInt.diagonal().real, 'b')
    plt.xlabel(dim_name[0])
    plt.ylabel("diag{Gamma}")

    # plt.savefig("0.90_Gamma.pdf")
    plt.show()

else:
    X = np.arange(0, Beta, Beta/N)
    Y = np.arange(0, Beta, Beta/N)
    X, Y = np.meshgrid(X, Y)

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    #for key in Quans:
    surf = ax.plot_surface(X, Y, GamInt.real, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
        #surf = ax.plot_surface(X, Y, GamMC[key][0].imag, rstride=1, cstride=1, cmap=cm.coolwarm,
            #linewidth=0, antialiased=False)
    # ax.set_zlim(-1.01, 1.01)

    ax.set_xlabel(dim_name[0])
    ax.set_ylabel(dim_name[1])
    ax.set_zlabel("Gamma")

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=Beta, aspect=5)

    plt.show()
