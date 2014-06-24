#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

fig = plt.figure()
ax = plt.subplot(111)

Quans=["SUMChi"]
key = "SUMChi"

SUMChi=[]
SUMChi.append(read_data.read_array("../../data/3D/bare_L4_0.50_1_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/3D/bare_L4_0.50_2_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/3D/bare_L4_0.50_3_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/3D/bare_L4_0.50_4_quantities.dat", Quans))

Order=np.arange(1, len(SUMChi)+1)
ChiL4=[]
for i in range(len(SUMChi)):
    ChiL4.append(sum(SUMChi[i][key][0]).real/N)

ax.plot(1.0/Order, ChiL4, marker='o', label="L=4, beta=0.50")

SUMChi=[]
SUMChi.append(read_data.read_array("../../data/3D/bare_L8_0.50_1_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/3D/bare_L8_0.50_2_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/3D/bare_L8_0.50_3_quantities.dat", Quans))

Order=np.arange(1, len(SUMChi)+1)
ChiL8=[]
for i in range( len(SUMChi)):
    ChiL8.append(sum(SUMChi[i][key][0]).real/N)

ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta=0.50")

ax.legend()

plt.xlabel("1/N")
plt.ylabel("Chi")

#plt.savefig("sum_chi_vs_tau_L8.pdf")
plt.savefig("Beta0.5_static_uniform_chi.pdf")
plt.show()
