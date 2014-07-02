#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data


def get_r(x,y,z,L):
    r = z*L**2+y*L+x
    return r

Beta = 0.80
L = 8

N = 64
omega = np.arange(0, N)
stag = (L**2+L+1)*(L/2)

Quans=["Denom"]
key="Denom"

Files=[]
Files.append(read_data.read_array("../../../../project/0.80_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)


for i in range(0, len(Files)):
    k=[]
    Denom = []
    for dx in range(0, L/2+1):
        for dy in range(0, L/2+1):
            for dz in range(0, L/2+1):
                k.append(np.sqrt(dx**2.0+dy**2.0+dz**2.0))
                Denom.append(Files[i][key][0][get_r(dx,dy,dz,L)][0].real)
    ax.plot(k, Denom, 'o', label="Beta="+str(Beta))

ax.legend()

plt.xlabel("|k|")
plt.ylabel("Denom(omega=0)")

plt.savefig("Beta0.8_L8_denominator_k.pdf")
plt.show()
