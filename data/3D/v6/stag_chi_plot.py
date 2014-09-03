#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.70
MxT = 128
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)

MnOrder = 3
MxOrder = 6
Quans=["ChiK0t"]
key="ChiK0t"

Files=[]
Files.append(read_data.read_array("L8_0.70_3/0.70_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.70_4/0.70_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.70_5/0.70_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.70_6/0.70_quantities.dat", Quans))

Order=np.arange(MnOrder, MxOrder+1)

L = 8
stag = (L**2+L+1)*(L/2)
#unif = 0

ChiL8=[]
for i in range(len(Files)):
    ChiL8.append(Files[i][key][0][stag].real)

ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta="+str(Beta))
########################################################################################

ax.legend()
ax.set_xlim((0.0, 1.0))

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.7_static_staggered_chi.pdf")
plt.show()
