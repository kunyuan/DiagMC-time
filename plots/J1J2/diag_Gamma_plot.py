#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.50
N = 64

tau = np.arange(0, Beta, Beta/N)

Quan2D=["Gamma1","Gamma2","Gamma3","Gamma4","Gamma5"]
#Quan2D=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Quan2D=["Gamma1","Gamma2","Gamma3"]
#Quan2D=["Gamma1","Gamma2"]
#Quan2D=["Gamma1"]
#Quan2D=[]

Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3","GammaDiag4","GammaDiag5"]
#Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3","GammaDiag4"]
#Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3"]
#Quan1D=["GammaDiag1","GammaDiag2"]
#Quan1D=["GammaDiag1"]
#Quan1D=[]
Quan=Quan2D+Quan1D

DiagGamma=[]
#DiagGamma.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_1_quantities.dat", Quan))
#DiagGamma.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_2_quantities.dat", Quan))
#DiagGamma.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_3_quantities.dat", Quan))
#DiagGamma.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_4_quantities.dat", Quan))
DiagGamma.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_5_quantities.dat", Quan))


fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(DiagGamma)):
    for key in Quan2D:
        ax.plot(tau, DiagGamma[i][key][0].diagonal().real, label=key)
        #ax.plot(tau, DiagGamma[i][key][0].diagonal().imag, label=key)
    for key in Quan1D:
        ax.plot(tau, DiagGamma[i][key][0].real, label=key)
        #ax.plot(tau, DiagGamma[i][key][0].imag, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

plt.savefig("diag_Gamma_5.pdf")
plt.show()