import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import read_data

is2d = True
#is2d = False

Beta = 0.5
N = 64

tau = np.arange(0, Beta, Beta/N)
#Gamt, dim_name = read_data.read_array("./0.90_Gam1.dat",scope="-1:")["Gamma"]
Gam0t, dim_name = read_data.read_array("./../0.90_Gam1.dat",scope="-1:")["Gamma"]
Gamt, dim_name = read_data.read_array("./../0.90_quantities.dat",["Gamma"])["Gamma"]

#tau = np.arange(0, Beta, Beta/N)
#SUMChit, dim_name = read_data.read_array("./1.00_quantities.dat",["SUMChi"])["SUMChi"]

#tau = np.arange(0, Beta, Beta/N)
#Denomt, dim_name = read_data.read_array("./1.00_quantities.dat",["Denom"])["Denom"]

tau = np.arange(0, Beta, Beta/N)
#Sigmat, dim_name = read_data.read_array("./bare_0.50_1_Sigma.dat",["Sigma"])["Sigma"]
#Sigmat, dim_name = read_data.read_array("./bare_0.50_2_Sigma.dat",["Sigma"])["Sigma"]
#Sigmat, dim_name = read_data.read_array("./1.00_quantities.dat",["Sigma"])["Sigma"]

if is2d==True:
    fig = plt.figure()

    plt.plot(tau, Gamt.diagonal().real, 'r', tau, Gam0t.diagonal().real,'b')
    plt.xlabel(dim_name[0])
    plt.ylabel("diag{Gamma}")

    #plt.plot(tau, SUMChit.real)
    #plt.xlabel(dim_name[0])
    #plt.ylabel("sum{Chi(r)}")

    #plt.plot(tau, Denomt[2,2,:].real)
    #plt.xlabel(dim_name[2])
    #plt.ylabel("Denom(0,0)")

    #plt.plot(tau, Sigmat.real)
    #plt.xlabel(dim_name[0])
    #plt.ylabel("Sigma")

    #plt.savefig("bare_Sigma_0.50_2.pdf")
    plt.savefig("Gamma_integral.pdf")
    
else:
    X = np.arange(0, Beta, Beta/N)
    Y = np.arange(0, Beta, Beta/N)
    X, Y = np.meshgrid(X, Y)
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_surface(
        X, Y, Gamt.real, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    # ax.set_zlim(-1.01, 1.01)

    ax.set_xlabel(dim_name[0])
    ax.set_ylabel(dim_name[1])
    ax.set_zlabel("Gamma")

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=Beta, aspect=5)

    plt.show()
