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

Quans=["StagChit"]
key = "StagChit"

Files=[]
Files.append(read_data.read_array("L8_0.90_1_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_2_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_3_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_4_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_5_quantities.dat", Quans))
Files.append(read_data.read_array("L8_0.90_6_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 8

fig = plt.figure() 
ax = plt.subplot(111)

for i in range(0, len(Files)):
    ax.plot(tau, Files[i][key][0].real, marker='o', label="L=8, Order "+str(i+1))
#ax.errorbar(0.05, 48.85, yerr=0.15, marker='*')
########################################################################################

ax.legend() 
plt.xlabel("tau")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.8_static_staggered_chi.pdf")
plt.show()
