#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

def get_r(x,y,z,L):
    r = z*L**2+y*L+x
    return r

Beta = 0.80
Pi = 3.1415926
N = 64
omega = np.arange(0, 2*Pi*N/Beta, 2*Pi/Beta)
L = 8
stag = (L**2+L+1)*(L/2)

Quans=["Denom"]
key="Denom"

Files=[]
Files.append(read_data.read_array("../../../../project/0.80_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)


for i in range(0, len(Files)):
    ax.plot(omega, Files[i][key][0][stag].real, 'o', label="Beta="+str(Beta))

ax.legend()

plt.xlabel("omega")
plt.ylabel("Denom(pi, pi, pi)")

plt.savefig("Beta0.8_L8_denominator.pdf")
plt.show()
