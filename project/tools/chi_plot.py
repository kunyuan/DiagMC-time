#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = np.arange(0, 4)

Quans=["Chi"]

Files=[]
Files.append(read_data.read_array("./0.50_0.10_4_quantities.dat", Quans))
Files.append(read_data.read_array("./0.50_0.50_4_quantities.dat", Quans))
Files.append(read_data.read_array("./0.50_1.00_4_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(Files)):
    for key in Quans:
        ax.plot(L, Files[i][key][0][0].real, label=key)

ax.legend()

plt.xlabel("dy")
plt.ylabel("Chi(0, dy)")

#plt.savefig("0.90_bold_Chi.pdf")
#plt.savefig("0.90_bare_Chi.pdf")
plt.show()
