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

Quans=["ChiKt0"]
key = "ChiKt0"

Files=[]
Files.append(read_data.read_array("L8_0.90_1_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_2_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_3_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_4_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_5_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_6_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 8

ChiL8=[]
for i in range(0, len(Files)):
    ChiL8.append(Files[i][key][0][L/2][L/2].real*L**2.0/3.0)
    print "L=8, Order ", i+1, ChiL8[i]

fig = plt.figure() 
ax = plt.subplot(111)


ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta="+str(Beta))
ax.errorbar(0.05, 48.85, yerr=0.15, marker='*')
########################################################################################

ax.legend() 
plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.8_static_staggered_chi.pdf")
plt.show()
