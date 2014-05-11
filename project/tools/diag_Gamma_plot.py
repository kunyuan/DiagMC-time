#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

#Quans=["Gamma"]
#Quans=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Quans=["Gamma1","Gamma2","Gamma3"]
#Quans=["Gamma1","Gamma"]

DiagGammas=[]
DiagGammas.append(read_data.read_array("../0.50_0.50_Gam1.dat", Quans))

Quan2D=["Gamma1","Gamma2","Gamma3","Gamma4","Gamma5"]
#Quan2D=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Quan2D=["Gamma1","Gamma2","Gamma3"]
#Quan2D=["Gamma1","Gamma2"]
#Quan2D=["Gamma1"]

Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3","GammaDiag4","GammaDiag5"]
#Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3","GammaDiag4"]
#Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3"]
#Quan1D=["GammaDiag1","GammaDiag2"]
#Quan1D=["GammaDiag1"]
#Quan1D=[]
Quan=Quan2D+Quan1D
#Quan=["Gamma"]
#print Quans

DiagGamma=[]
DiagGamma.append(read_data.read_array("./0.50_4_bold_quantities.dat", Quan))
#DiagGamma.append(read_data.read_array("./0.90_4_bare_quantities.dat", Quan))


fig = plt.figure()
ax = plt.subplot(111)

#for i in range(len(DiagGammas)):
    #for key in Quans:
        #ax.plot(tau, DiagGammas[i][key][0].diagonal().real, label=key)
        ##ax.plot(tau, DiagGamma[i][key][0].diagonal().imag, label=key)

for i in range(len(DiagGamma)):
    for key in Quan2D:
        ax.plot(tau, ((DiagGamma[i][key][0].diagonal()+DiagGamma[i][key][0].diagonal()[::-1])/2).real, label=key)
        #ax.plot(tau, DiagGamma[i][key][0].diagonal().imag, label=key)
    for key in Quan1D:
        ax.plot(tau, ((DiagGamma[i][key][0]+DiagGamma[i][key][0][::-1])/2).real, label=key)
        #ax.plot(tau, DiagGamma[i][key][0].real, label=key)
        #ax.plot(tau, DiagGamma[i][key][0].imag, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
