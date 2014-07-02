#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.80

MxT = 64
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)


Quans=["ChiK"]
key="ChiK"


########################### L=8 ########################################################
Files=[]
Files.append(read_data.read_array("../../../../project/0.80_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 8
stag = (L**2+L+1)*(L/2)

ChiL8=[]
for i in range(0, len(Files)):
    ChiL8.append(Files[i][key][0][stag].real)

ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta="+str(Beta))
########################################################################################

ax.legend()
ax.set_xlim(-0.05, 1.05)

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.8_static_staggered_chi.pdf")
plt.show()
