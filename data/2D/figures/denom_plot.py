#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = 8
Beta = 0.90
N = 64
tau = np.arange(0, Beta, Beta/N)

Quans=["Denom"]

Files=[]
#Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L4_0.90_1_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L4_0.90_2_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L4_0.90_3_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L4_0.90_4_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L4_0.90_5_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L4_0.90_6_quantities.dat", Quans))

Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L8_0.90_1_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L8_0.90_2_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L8_0.90_3_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L8_0.90_4_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L8_0.90_5_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/2D/cmp_L_dependence/L8_0.90_6_quantities.dat", Quans))




fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(Files)):
    for key in Quans:
        ax.plot(tau, Files[i][key][0][L/2][L/2].real, marker='o', label="Order "+str(i+1))

ax.legend()

plt.xlabel("omega")
plt.ylabel("Denom(pi, pi)")

#plt.savefig("Beta0.9_denominator_L4.pdf")
plt.savefig("Beta0.9_denominator_L8.pdf")

plt.show()
