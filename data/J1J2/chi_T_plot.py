#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

MxT = 128

fig = plt.figure()
ax = plt.subplot(111)

Quans=["ChiK"]
key = "ChiK"

Files5=[]
Files5.append(read_data.read_array("0.67/0.67_quantities.dat", Quans))
Files5.append(read_data.read_array("1.00/1.00_quantities.dat", Quans))
Files5.append(read_data.read_array("1.50/1.50_quantities.dat", Quans))
Files5.append(read_data.read_array("1.75/1.75_quantities.dat", Quans))
Files5.append(read_data.read_array("2.00/2.00_quantities.dat", Quans))
Files5.append(read_data.read_array("2.25/2.25_quantities.dat", Quans))

Files0=[]
Files0.append(read_data.read_array("1.00_mf/1.00_quantities.dat", Quans))
Files0.append(read_data.read_array("1.50_mf/1.50_quantities.dat", Quans))
Files0.append(read_data.read_array("2.00_mf/2.00_quantities.dat", Quans))
Files0.append(read_data.read_array("2.50_mf/2.50_quantities.dat", Quans))
Files0.append(read_data.read_array("2.75_mf/2.75_quantities.dat", Quans))

T5 =[]
T5.append(1.0/0.67)
T5.append(1.0/1.00)
T5.append(1.0/1.50)
T5.append(1.0/1.75)
T5.append(1.0/2.00)
T5.append(1.0/2.25)

T0 =[]
T0.append(1.0/1.00)
T0.append(1.0/1.50)
T0.append(1.0/2.00)
T0.append(1.0/2.50)
T0.append(1.0/2.75)

L = 16
stag = (L+1)*(L/2)
unif = 0

ChiL8U=[]
ChiL8S=[]
for i in range(len(Files5)):
    ChiL8U.append(Files5[i][key][0][unif].real/T5[i])
    ChiL8S.append(Files5[i][key][0][stag].real/T5[i])

ax.plot(T5, ChiL8U, marker='o', label="uniform, L=16, J2/J1=0.5, Order= 5")
#ax.plot(T5, ChiL8S, marker='o', label="staggered, L=16, J2/J1=0.5, Order=5")

ChiL8U=[]
ChiL8S=[]
for i in range(len(Files0)):
    ChiL8U.append(Files0[i][key][0][unif].real/T0[i])
    ChiL8S.append(Files0[i][key][0][stag].real/T0[i])

ax.plot(T0, ChiL8U, marker='o', label="uniform, L=16, J2/J1=0.5, Order=0")
#ax.plot(T0, ChiL8S, marker='o', label="staggered, L=16, J2/J1=0.5, Order=0")


########################################################################################

ax.legend()
#ax.set_xlim(0.0)
#ax.set_ylim(0.0)

plt.xlabel(r"T")
plt.ylabel(r"$\chi_u$")
#plt.ylabel(r"$\chi_s$")

plt.savefig("static_uniform_chi_T.pdf")
#plt.savefig("static_staggered_chi_T.pdf")
#plt.savefig("static_chi_T.pdf")
plt.show()
