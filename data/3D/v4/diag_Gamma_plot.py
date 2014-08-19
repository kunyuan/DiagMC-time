#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.70
N = 128

tau = np.arange(0, Beta, Beta/N)

Quan2D=["Gamma1","Gamma2","Gamma3","Gamma4","Gamma5","Gamma6"]
#Quan2D=["Gamma1","Gamma2","Gamma3","Gamma4","Gamma5"]
#Quan2D=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Quan2D=["Gamma3","Gamma"]
#Quan2D=["Gamma1","Gamma2"]
#Quan2D=["Gamma1"]
#Quan2D=["Gamma"]
#Quan2D=[]

Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3","GammaDiag4","GammaDiag5","GammaDiag6"]
#Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3","GammaDiag4","GammaDiag5"]
#Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3","GammaDiag4"]
#Quan1D=["GammaDiag1","GammaDiag2","GammaDiag3"]
#Quan1D=["GammaDiag1","GammaDiag2"]
#Quan1D=["GammaDiag1"]
#Quan1D=[]
Quan=Quan2D+Quan1D

DiagGamma=[]
DiagGamma.append(read_data.read_array("L8_0.70_6/0.70_quantities.dat", Quan))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(DiagGamma)):
    for key in Quan2D:
        ax.plot(tau, DiagGamma[i][key][0].diagonal().real, marker='o', label=key)
        #ax.plot(tau, DiagGamma[i][key][0].diagonal().imag, marker='o', label=key)
    for key in Quan1D:
        ax.plot(tau, DiagGamma[i][key][0].real, marker='*', label=key)
        ##ax.plot(tau, DiagGamma[i][key][0].imag, marker='*', label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

plt.savefig("Beta0.70_6_diag_Gamma.pdf")
plt.show()
