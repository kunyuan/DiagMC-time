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
Order=np.arange(1, 7)

Quans=["SUMChi"]

SUMChi=[]
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L4_0.90_1_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L4_0.90_2_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L4_0.90_3_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L4_0.90_4_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L4_0.90_5_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L4_0.90_6_quantities.dat", Quans))

key = "SUMChi"
ChiL4=[]
for i in range(len(SUMChi)):
    ChiL4.append(sum(SUMChi[i][key][0]).real/N)

SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L8_0.90_1_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L8_0.90_2_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L8_0.90_3_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L8_0.90_4_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L8_0.90_5_quantities.dat", Quans))
SUMChi.append(read_data.read_array("../../data/cmp_L_dependence/L8_0.90_6_quantities.dat", Quans))

ChiL8=[]
for i in range(len(ChiL4), len(SUMChi)):
    ChiL8.append(sum(SUMChi[i][key][0]).real/N)

fig = plt.figure()
ax = plt.subplot(111)

#for i in range(len(SUMChi)):
    #for key in Quans:
        #ax.plot(tau, SUMChi[i][key][0].real, label="Order "+str(i+1))

ax.plot(1.0/Order, ChiL4, marker='o', label="L=4, beta=0.90")
ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta=0.90")

ax.legend()

plt.xlabel("1/N")
plt.ylabel("Chi")

#plt.savefig("sum_chi_vs_tau_L8.pdf")
#plt.savefig("static_uniform_chi.pdf")
plt.show()
