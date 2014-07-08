#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.70
N = 64
tau = np.arange(0, Beta, Beta/N)

L = 16
stag = (L+1)*(L/2)

Quans=["Denom"]

Files=[]
Files.append(read_data.read_array("./L16_0.70_1_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_2_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_3_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_4_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_5_bare_quantities.dat", Quans))


fig = plt.figure()
ax = plt.subplot(111)


for i in range(len(Files)):
    for key in Quans:
        ax.plot(tau, Files[i][key][0][stag].real, marker='o', label="Order "+str(i+1))

ax.legend()

plt.xlabel("omega")
plt.ylabel("Denom(pi, pi)")

plt.savefig("Beta0.7_L16_denominator.pdf")
plt.show()
