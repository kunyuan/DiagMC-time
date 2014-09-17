#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

def get_r(x,y,z,L):
    r = z*L**2+y*L+x
    return r

Quans=["GammaR"]

##################################L=4#################################
L = 8
Lx = np.arange(0, L)
Ly = np.arange(0, L)
Lz = np.arange(0, L)

r=[]
for x in Lx:
    for y in Ly:
        for z in Lz:
            r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

Files=[]
Files.append(read_data.read_array("project/0.50_quantities.dat", Quans))
Files.append(read_data.read_array("v7.3.1/0.50_quantities.dat", Quans))

fig = plt.figure() 
ax = plt.subplot(111)

for key in Quans:
    for i in range(len(Files)):
        ax.plot(r, Files[i][key][0].real, 'o', label=key)

ax.legend()

plt.xlabel("dr") 
plt.ylabel("Gamma(dr)")

plt.savefig("Beta0.5_L8_Gamma_r.pdf")
plt.show()

