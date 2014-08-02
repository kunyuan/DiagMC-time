#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.50
N = 128

tau = np.arange(0, Beta, Beta/N)
target=["Sigma2", "Sigma3", "Sigma"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("./bare_L4_0.50_3_quantities.dat",target))
BoldSigma.append(read_data.read_array("./Sigma_Order3_Nikolay.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for key in target:
    ax.plot(tau, BoldSigma[0][key][0].real, marker='^', label="{0}, Yuan's".format(key))

for key in target:
    ax.plot(tau, -1.0*BoldSigma[1][key][0].real, marker='o', label="{0}, Nikolay's".format(key))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
