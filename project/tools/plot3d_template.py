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

#Gamma = []
#G = read_data.read_array("./../1.00_Gam1.dat")
#Gamma.append(G)

#G = read_data.read_array("./../1_1.00_1_Gam_typ.dat")
#Gamma.append(G)
#G = read_data.read_array("./../2_1.00_1_Gam_typ.dat")
#Gamma.append(G)
#G = read_data.read_array("./../3_1.00_1_Gam_typ.dat")
#Gamma.append(G)
#G = read_data.read_array("./../4_1.00_1_Gam_typ.dat")
#Gamma.append(G)

#filename = ["Gamma1", "Gamma2", "Gamma3", "Gamma4", "Gamma5"]

#GamMC = []
#GamInt = []
#for i in range(5):
    #G0, dim_name = Gamma[0][filename[i]]
    #GamInt.append(G0)

    #Gamt=[]
    #Gt, dim_name = Gamma[1][filename[i]]
    #Gamt.append(Gt)
    #Gt, dim_name = Gamma[2][filename[i]]
    #Gamt.append(Gt)
    #Gt, dim_name = Gamma[3][filename[i]]
    #Gamt.append(Gt)
    #Gt, dim_name = Gamma[4][filename[i]]
    #Gamt.append(Gt)

    #mGamt = (Gamt[0]+Gamt[1]+Gamt[2]+Gamt[3])/4.0
    #GamMC.append(mGamt)

#tau = np.arange(0, Beta, Beta/N)
#SUMChit, dim_name = read_data.read_array("./1.00_quantities.dat",["SUMChi"])["SUMChi"]

#tau = np.arange(0, Beta, Beta/N)
#Denomt, dim_name = read_data.read_array("./1.00_quantities.dat",["Denom"])["Denom"]

tau = np.arange(0, Beta, Beta/N)
#Sigmat, dim_name = read_data.read_array("./bare_0.50_1_Sigma.dat",["Sigma"])["Sigma"]
#Sigmat, dim_name = read_data.read_array("./bare_0.50_2_Sigma.dat",["Sigma"])["Sigma"]
#Sigmat, dim_name = read_data.read_array("./1.00_quantities.dat",["Sigma"])["Sigma"]

if is2d is True:
    fig = plt.figure()
    plt.plot(tau, GamMC.diagonal().real, 'r', 
            #tau, GamMC[1].diagonal().real, 'r', 
            #tau, GamMC[2].diagonal().real,'b', 
            #tau, GamMC[3].diagonal().real, 'b*',
            #tau, GamMC[4].diagonal().real,'g--', 
            tau, GamInt.diagonal().real, 'b')
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
    # plt.savefig("Gamma_integral.pdf")
    plt.show()

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
