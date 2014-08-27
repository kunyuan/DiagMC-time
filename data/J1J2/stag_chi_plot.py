#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 1.50
MxT = 128
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)

MxOrder = 5
Quans=["ChiK0t1", "ChiK0t2", "ChiK0t3", "ChiK0t4", "ChiK0t5"]

Files=[]
Files.append(read_data.read_array("1.50/1.50_order_quantities.dat", Quans))

Order=np.arange(1, MxOrder+1)

L = 16
stag = (L+1)*(L/2)
unif = 0

ChiL8U=[]
ChiL8S=[]
for key in Quans:
    ChiL8U.append(Files[0][key][0][unif].real)
    ChiL8S.append(Files[0][key][0][stag].real)
#print ChiL8[MxOrder-1]

ax.plot(1.0/Order, ChiL8U, marker='o', label="uniform, L=16, J2/J1=0.5, beta="+str(Beta))
ax.plot(1.0/Order, ChiL8S, marker='o', label="staggered, L=16, J2/J1=0.5, beta="+str(Beta))
########################################################################################

ax.legend()
ax.set_xlim((0.0, 1.0))

plt.xlabel("1/N")
plt.ylabel("susceptibility")

plt.savefig("Beta1.5_static_chi.pdf")
plt.show()
