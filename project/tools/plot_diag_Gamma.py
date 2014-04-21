import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import read_data

is2d = True
#is2d = False

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)
GamInt, dim_name = read_data.read_array("./../0.90_Gam1.dat")["Gamma"]

GamMC, dim_name = read_data.read_array("./../0.90_quantities.dat")["Gamma"]
#Gamma.append(GamMC)
#Gamma=[]
#GamMC, dim_name = read_data.read_array("./../data/bold_0.90_2_Gammas.dat")["Gamma2"]
#Gamma.append(GamMC)
#GamMC, dim_name = read_data.read_array("./../data/bold_0.90_3_Gammas.dat")["Gamma3"]
#Gamma.append(GamMC)
#GamMC, dim_name = read_data.read_array("./../data/bold_0.90_4_Gammas.dat")["Gamma4"]
#Gamma.append(GamMC)


if is2d is True:
    fig = plt.figure()
    plt.plot(tau, GamMC.diagonal().real, 'r', 
            tau, GamInt.diagonal().real, 'b')
    #plt.plot(
            #tau, Gamma[0].diagonal().real, 'r',
            #tau, Gamma[1].diagonal().real, 'b',
            #tau, Gamma[2].diagonal().real, 'g',
            #tau, Gamma[3].diagonal().real, 'r--',
            #)

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

    surf = ax.plot_surface(
        X, Y, GamMC.real, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    # ax.set_zlim(-1.01, 1.01)

    ax.set_xlabel(dim_name[0])
    ax.set_ylabel(dim_name[1])
    ax.set_zlabel("Gamma")

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=Beta, aspect=5)

    plt.show()
