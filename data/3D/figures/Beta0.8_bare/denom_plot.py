#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.80
N = 64
tau = np.arange(0, Beta, Beta/N)

Quans=["Denom"]

Files=[]
Files.append(read_data.read_array("./0.80_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.80_d_quantities.dat", Quans))


fig = plt.figure()
ax = plt.subplot(111)

L = 8
stag = (L**2+L+1)*(L/2)


for i in range(len(Files)):
    for key in Quans:
        ax.plot(tau, Files[i][key][0][stag].real, marker='o', label="Order "+str(i+1))

ax.legend()

plt.xlabel("omega")
plt.ylabel("Denom(pi, pi)")

plt.savefig("Beta0.8_L4_denominator.pdf")
#plt.savefig("denominator_L8.pdf")
plt.show()
