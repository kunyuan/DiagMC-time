#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.50
N = 64

Quans=["Sigma"]

tau = np.arange(0, Beta, Beta/N)

BoldSigma=[]
BoldSigma.append(read_data.read_array("../0.50_0.40_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in Quans:
        ax.plot(tau, BoldSigma[i][key][0].real, label="{0}, Order{1}".format(key, i+2))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
