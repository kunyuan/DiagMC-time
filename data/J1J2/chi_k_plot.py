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
Files.append(read_data.read_array("1.50_quantities.dat", Quans))

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


    ax.plot(path, realChiK, marker='o', label=r"$\Gamma-X-M-\Gamma$ Order "+str(i+3))


ax.legend()

plt.ylabel("Chi(k)")
plt.savefig("Beta1.50_L16_Chi_k.pdf")

plt.show()

