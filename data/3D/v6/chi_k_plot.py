#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

def get_r(x,y,z,L):
    r = z*L**2+y*L+x
    return r

N = 128.0

Quans=["ChiK"]
key="ChiK"

fig = plt.figure()
ax = plt.subplot(111)


L = 8
Lx = np.arange(0, L/2)
Ly = np.arange(0, L/2)
Lz = np.arange(0, L/2)

Order =5

Beta=[]
Beta.append(0.5)
Beta.append(0.6)
Beta.append(0.7)

Files=[]
Files.append(read_data.read_array("L8_0.50_"+str(Order)+"/0.50_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.60_"+str(Order)+"/0.60_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.70_"+str(Order)+"/0.70_quantities.dat", Quans))


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
        r = get_r(0,y,0,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath+1, -0.22500, r"$X$")

#######X-M
    for x in Lx:
        ipath = ipath + 1
        path.append(ipath)
        r = get_r(x,L/2,0,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath+1, -0.2500, r"$M$")


#######M-Gamma
    for x in Lx:
        ipath = ipath + 1
        path.append(ipath)
        y = x
        r = get_r(L/2-x,L/2-y,0,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath+1, -0.2500, r"$\Gamma$")


#######Gamma-R
    for x in Lx:
        ipath = ipath + 1
        path.append(ipath)
        y = x
        z = x
        r = get_r(x,y,z,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)
    ax.text(ipath+1, -0.2500, r"$R$")


#######R-X
    for x in Lx:
        ipath = ipath + 1
        path.append(ipath)
        z = x
        r = get_r(L/2-x,L/2,L/2-z,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath, -0.2500, r"$X$")

    ax.plot(path, realChiK, marker='o', label=r"$\Gamma-X-M-\Gamma-R-X$ Beta "+str(Beta[i]))

#######M-R
    path = []
    realChiK = []
    imagChiK = []
    ipath = ipath + 1
    ax.text(ipath+1, -0.2500, r"$M$")
    for z in Lz:
        ipath = ipath + 1
        path.append(ipath)
        r = get_r(L/2,L/2,z,L)
        realChiK.append(Files[i][key][0][r].real)
        imagChiK.append(Files[i][key][0][r].imag)

    ipath = ipath + 1
    path.append(ipath)
    r = get_r(L/2,L/2,L/2,L)
    realChiK.append(Files[i][key][0][r].real)
    imagChiK.append(Files[i][key][0][r].imag)

    ax.text(ipath, -0.2500, r"$R$")
    ax.plot(path, realChiK, marker='*', label=r"$M-R$ "+str(Beta[i]))


ax.legend()

plt.ylabel("Chi(k)")

plt.savefig("L8_0.50_0.70_Chi_k.pdf")

plt.show()

