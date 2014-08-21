#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

def get_r(x,y,L):
    r = y*L+x
    return r

N = 128

Quans=["ChiK"]
key="ChiK"

fig = plt.figure()
ax = plt.subplot(111)


L = 16
Lx = np.arange(0, L/2)
Ly = np.arange(0, L/2)

Files=[]
#Files.append(read_data.read_array("1.50/1.50_quantities.dat", Quans))
#Files.append(read_data.read_array("1.75/1.75_quantities.dat", Quans))
#Files.append(read_data.read_array("2.00/2.00_quantities.dat", Quans))
Files.append(read_data.read_array("1.00_mf/1.00_quantities.dat", Quans))
Files.append(read_data.read_array("1.50_mf/1.50_quantities.dat", Quans))
Files.append(read_data.read_array("2.00_mf/2.00_quantities.dat", Quans))
Files.append(read_data.read_array("2.50_mf/2.50_quantities.dat", Quans))
Files.append(read_data.read_array("2.75_mf/2.75_quantities.dat", Quans))

Beta=[]
Beta.append(1.00)
Beta.append(1.50)
Beta.append(2.00)
Beta.append(2.50)
Beta.append(2.75)

for i in range(len(Files)):

#######Gamma-X
    path = []
    realChiK = []
    imagChiK = []
    ipath = 0
    ax.text(ipath+1, -0.225000, r"$\Gamma$")
    for y in Ly:
        ipath = ipath + 1
        path.append(ipath)
        r = get_r(0,y,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath+1, -0.22500, r"$X$")

#######X-M
    for x in Lx:
        ipath = ipath + 1
        path.append(ipath)
        r = get_r(x,L/2,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath+1, -0.2500, r"$M$")


#######M-Gamma
    for x in Lx:
        ipath = ipath + 1
        path.append(ipath)
        y = x
        r = get_r(L/2-x,L/2-y,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath+1, -0.2500, r"$\Gamma$")


    ax.plot(path, realChiK, marker='o', label=r"$\Gamma-X-M-\Gamma$ Beta="+str(Beta[i]))


ax.legend()

plt.ylabel("Chi(k)")
plt.savefig("L16_J1J2_Chi_k.pdf")

plt.show()

