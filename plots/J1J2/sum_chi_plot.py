#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.50
N = 64

tau = np.arange(0, Beta, Beta/N)
Order=np.arange(1, 6)

Quans=["SUMChi"]
key = "SUMChi"

SUMChi=[]
SUMChi.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_1_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_2_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_3_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_4_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_5_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

#ChiL4=[]
#for i in range(len(SUMChi)):
    #ChiL4.append(sum(SUMChi[i][key][0]).real/N)
#ax.plot(1.0/Order, ChiL4, marker='o', label="L=4, beta=0.90")

for key in Quans:
    for i in range(len(SUMChi)):
        ax.plot(tau, SUMChi[i][key][0].real, marker='*', label='L=4 Order '+str(i+1))

ax.legend()

#plt.xlabel("1/N")
#plt.ylabel("Chi")
#plt.savefig("static_uniform_chi.pdf")

plt.xlabel("tau")
plt.ylabel("Chi")
plt.savefig("chi_tau.pdf")

plt.show()
