#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = 8
Beta = 0.70
strBeta="0.70"
N = 64
tau = np.arange(0, Beta, Beta/N)

Quans=["Denom"]

Files=[]
Files.append(read_data.read_array("./bare_L"+str(L)+"_"+strBeta+"_1_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L"+str(L)+"_"+strBeta+"_2_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L"+str(L)+"_"+strBeta+"_3_quantities.dat", Quans))


fig = plt.figure()
ax = plt.subplot(111)

stag = (L**2+L+1)*(L/2)


for i in range(len(Files)):
    for key in Quans:
        ax.plot(tau, Files[i][key][0][stag].real, marker='o', label="Order "+str(i+1))
        ax.plot(tau, Files[i][key][0][stag].imag, marker='o', label="Order "+str(i+1))
        print Files[i][key][0][stag][0]

ax.legend()

plt.xlabel("omega")
plt.ylabel("Denom(pi, pi)")

plt.savefig("Beta"+str(Beta)+"_L"+str(L)+"_denominator.pdf")
plt.show()
