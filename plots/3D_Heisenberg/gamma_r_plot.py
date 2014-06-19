#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data


L = 4
Lx = np.arange(0, L)
Ly = np.arange(0, L)
Lz = np.arange(0, L)

r=[]
for x in Lx:
    for y in Ly:
        for z in Lz:
            r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

x = 0.0

N = 64.0

Quans=["GammaR"]

Files=[]
Files.append(read_data.read_array("../../project/0.50_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)


for key in Quans:
    for i in range(len(Files)):
        ax.plot(r, Files[i][key][0].real, marker='o', label=key+" Order"+str(i+1))

ax.legend()


plt.xlabel("dr") 
plt.ylabel("Gamma(dr)")

plt.show()
