#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

#Beta = 0.50
Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)
#target=["G0","G1","G2"]
target=["G"]

BoldSigma=[]
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_1_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_2_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_3_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_4_quantities.dat",target))
BoldSigma.append(read_data.read_array("../../project/0.80_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        ax.plot(tau, BoldSigma[i][key][0].real, label="{0}, Order{1}".format(key, i+1))

ax.legend()

plt.xlabel("tau")
plt.ylabel("G")

#plt.savefig("Beta0.5_L4_G.pdf")
plt.savefig("Beta0.9_L4_G.pdf")
plt.show()
