#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

#Gamma=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Gamma=["Gamma1","Gamma2","Gamma3"]
#Gamma=["Gamma1","Gamma2","Gamma"]
#Gamma=["Gamma"]



DiagGamma=[]
DiagGamma.append(read_data.read_array("./0.90_1_bold_quantities.dat", Quans))
DiagGamma.append(read_data.read_array("./0.90_2_bold_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(DiagGamma)):
    for key in Quans:
        ax.plot(tau, DiagGamma[i][key][0].diagonal().real, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
