import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

is2d = True
#is2d = False

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

BareChi=[]
Chi,dim_name=read_data.read_array("./bare_0.90_2_Chi_sum.dat")["SUMChi"]
BareChi.append(Chi)
Chi,dim_name=read_data.read_array("./bare_0.90_3_Chi_sum.dat")["SUMChi"]
BareChi.append(Chi)
Chi,dim_name=read_data.read_array("./bare_0.90_4_Chi_sum.dat")["SUMChi"]
BareChi.append(Chi)

BoldChi=[]
Chi,dim_name=read_data.read_array("../data/bold_0.90_2_quantities.dat")["SUMChi"]
BoldChi.append(Chi)
Chi,dim_name=read_data.read_array("../data/bold_0.90_3_quantities.dat")["SUMChi"]
BoldChi.append(Chi)
Chi,dim_name=read_data.read_array("../data/bold_0.90_4_quantities.dat")["SUMChi"]
BoldChi.append(Chi)

if is2d is True:
    fig = plt.figure()
    plt.plot(
            tau, BareChi[0].real, 'r',
            tau, BoldChi[0].real, 'b',
            tau, BareChi[1].real, 'r*',
            tau, BoldChi[1].real, 'b*',
            tau, BareChi[2].real, 'ro',
            tau, BoldChi[2].real, 'bo')
    plt.xlabel(dim_name[0])
    plt.ylabel("sum{Chi(r)}")

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
