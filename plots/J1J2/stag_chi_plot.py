#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

#Beta = 0.50
Beta = 0.70

MxT = 64
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)


Quans=["ChiK"]
key = "ChiK"

########################### L=16 ########################################################
Files=[]
Files.append(read_data.read_array("../../data/J1J2/L16_0.70_1_bare_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/J1J2/L16_0.70_2_bare_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/J1J2/L16_0.70_3_bare_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/J1J2/L16_0.70_4_bare_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 16
stag = (L+1)*(L/2)

ChiL16=[]
for i in range(0, len(Files)):
    ChiL16.append(Files[i][key][0][stag].real)

ax.plot(1.0/Order, ChiL16, marker='o', label="L="+str(L)+", beta="+str(Beta))
########################################################################################

ax.legend()
ax.set_xlim(-0.05, 1.05)

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.7_static_staggered_chi.pdf")
plt.show()
