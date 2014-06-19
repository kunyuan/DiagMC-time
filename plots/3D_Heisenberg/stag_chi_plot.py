#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.50
MxT = 64

tau = np.arange(0, Beta, Beta/N)
Order=np.arange(1, 7)

L = 4

stag = (L+1)*(L/2)

Quans=["ChiK"]

Files=[]
Files.append(read_data.read_array("../../data/3D/cmp_highT/L4_0.50_1_quantitites.dat", Quans))
Files.append(read_data.read_array("../../data/3D/cmp_highT/L4_0.50_2_quantitites.dat", Quans))
Files.append(read_data.read_array("../../data/3D/cmp_highT/L4_0.50_3_quantitites.dat", Quans))
Files.append(read_data.read_array("../../data/3D/cmp_highT/L4_0.50_4_quantitites.dat", Quans))

Files.append(read_data.read_array("../../data/3D/cmp_highT/L8_0.50_1_quantitites.dat", Quans))
Files.append(read_data.read_array("../../data/3D/cmp_highT/L8_0.50_2_quantitites.dat", Quans))
Files.append(read_data.read_array("../../data/3D/cmp_highT/L8_0.50_3_quantitites.dat", Quans))
Files.append(read_data.read_array("../../data/3D/cmp_highT/L8_0.50_4_quantitites.dat", Quans))

key = "ChiK"
ChiL4=[]
for i in range(len(Files)):
    ChiL4.append(Files[i][key][stag].real)

fig = plt.figure()
ax = plt.subplot(111)

ax.plot(1.0/Order, ChiL4, marker='o', label="L=4, beta=0.50")
ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta=0.50")

ax.legend()

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

#plt.savefig("sum_chi_vs_tau_L8.pdf")
#plt.savefig("static_uniform_chi.pdf")
plt.show()
