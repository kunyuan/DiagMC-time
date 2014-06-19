#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = np.arange(0, 4)
Beta = 0.90
N = 64
tau = np.arange(0, Beta, Beta/N)

Quans=["Denom"]

Files=[]
#Files.append(read_data.read_array("./L8_0.50_3_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L4_0.90_1_quantities.dat", Quans))
Files.append(read_data.read_array("./L4_0.90_2_quantities.dat", Quans))
Files.append(read_data.read_array("./L4_0.90_3_quantities.dat", Quans))
Files.append(read_data.read_array("./L4_0.90_4_quantities.dat", Quans))
Files.append(read_data.read_array("./L4_0.90_5_quantities.dat", Quans))

#Files.append(read_data.read_array("./L8_0.90_1_quantities.dat", Quans))
#Files.append(read_data.read_array("./L8_0.90_2_quantities.dat", Quans))
#Files.append(read_data.read_array("./L8_0.90_3_quantities.dat", Quans))
#Files.append(read_data.read_array("./L8_0.90_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./L8_0.90_5_quantities.dat", Quans))

#Files.append(read_data.read_array("./0.90_4_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_5_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_4_bare_quantities.dat", Quans))

#Files.append(read_data.read_array("./0.50_0.10_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.50_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.50_1.00_4_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(Files)):
    for key in Quans:
        ax.plot(tau, Files[i][key][0][2][2].real, marker='o', label="Order "+str(i+1))

ax.legend()

plt.xlabel("omega")
plt.ylabel("Denom(pi, pi)")

plt.savefig("denominator_L4.pdf")
#plt.savefig("denominator_L8.pdf")
#plt.show()
