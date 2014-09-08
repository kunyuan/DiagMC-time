#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 2.00
MxT = 128
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)

MnOrder = 1
MxOrder = 5
Quans=["ChiK0t"]
key="ChiK0t"

Files=[]
Files.append(read_data.read_array("L16_2.00_0/2.00_quantities.dat", Quans))
Files.append(read_data.read_array("L16_2.00_2/2.00_quantities.dat", Quans))
Files.append(read_data.read_array("L16_2.00_3/2.00_quantities.dat", Quans))
Files.append(read_data.read_array("L16_2.00_4/2.00_quantities.dat", Quans))
Files.append(read_data.read_array("L16_2.00_5/2.00_quantities.dat", Quans))

Order=np.arange(MnOrder, MxOrder+1)

L = 16 
#stag = (L+1)*(L/2)
unif = 0

Chi=[]
for i in range(len(Files)):
    Chi.append(Files[i][key][0][unif].real)

ax.plot(1.0/Order, Chi, marker='o', label="L=16, beta="+str(Beta))
########################################################################################

ax.legend()
ax.set_xlim((0.0, 1.0))

plt.xlabel("1/N")
plt.ylabel("uniform susceptibility")

plt.savefig("Beta2.0_static_uniform_chi.pdf")
#plt.savefig("Beta1.5_static_uniform_chi.pdf")
plt.show()
