#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = np.arange(0, 16)
Beta = 0.90
N = 64
tau = np.arange(0, Beta, Beta/N)

Quans=["Denom"]

Files=[]
Files.append(read_data.read_array("../0.90_quantities.dat", Quans))
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
        for j in range(0, 16):
            for k in range(0, 16):
                ax.plot(tau, Files[i][key][0][j][k].real, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("Denom(pi/2, pi/2)")

plt.show()
