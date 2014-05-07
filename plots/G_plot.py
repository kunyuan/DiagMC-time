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
#target=["G0","G1","G2"]
target=["G"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("../0.90_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        ax.plot(tau, BoldSigma[i][key][0].real, label="{0}, Order{1}".format(key, i))
        ax.plot(tau, BoldSigma[i][key][0][::-1].real, label="{0}, Order{1}".format(key, i))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
