#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.50
MxT = 64
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)


Quans=["ChiK"]

Files=[]
Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_1_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_2_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_3_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_4_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 4
stag = (L**2+L+1)*(L/2)

key = "ChiK"
ChiL4=[]
for i in range(0, len(Files)):
    ChiL4.append(Files[i][key][0][stag].real)
ax.plot(1.0/Order, ChiL4, marker='o', label="L=4, beta=0.50")

Files=[]
Files.append(read_data.read_array("../../data/3D/bare_L8_0.50_1_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L8_0.50_2_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L8_0.50_3_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 8
stag = (L**2+L+1)*(L/2)

ChiL8=[]
for i in range(0, len(Files)):
    ChiL8.append(Files[i][key][0][stag].real)
ax.plot(1.0/Order, ChiL8, marker='^', label="L=8, beta=0.50")

ax.plot(0.001, 1.99598, marker='*', label="beta=0.50, high-T expansion")

ax.legend()

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.5_static_staggered_chi.pdf")
plt.show()
