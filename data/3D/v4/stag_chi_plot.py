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

MxOrder = 6
Quans=["ChiK0t1", "ChiK0t2", "ChiK0t3", "ChiK0t4", "ChiK0t5", "ChiK0t6"]
#Quans=["ChiK0t1", "ChiK0t2", "ChiK0t3", "ChiK0t4", "ChiK0t5"]

Files=[]
Files.append(read_data.read_array("L8_0.70_6/0.70_order_quantities.dat", Quans))
#Files.append(read_data.read_array("L8_0.60_5/0.60_order_quantities.dat", Quans))
#Files.append(read_data.read_array("L8_0.50_5/0.50_order_quantities.dat", Quans))

Order=np.arange(1, MxOrder+1)

L = 8
stag = (L**2+L+1)*(L/2)
#unif = 0

ChiL8=[]
for key in Quans:
    ChiL8.append(Files[0][key][0][stag].real)
print ChiL8[MxOrder-1]

ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta="+str(Beta))
########################################################################################

ax.legend()
ax.set_xlim((0.0, 1.0))

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.7_static_staggered_chi.pdf")
plt.show()
