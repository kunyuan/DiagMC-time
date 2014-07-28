#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.80
N = 64
L = 4

tau = np.arange(0, Beta, Beta/N)
target=["W"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("./0.80_quantities.dat",target))
BoldSigma.append(read_data.read_array("./0.80_normal_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        ax.plot(tau, BoldSigma[i][key][0][0,:].real, label="{0}, Order{1}".format(key, (0,0)))
        ax.plot(tau, BoldSigma[i][key][0][1,:].real, label="{0}, Order{1}".format(key, (1,0)))

ax.legend()

plt.xlabel("tau")
plt.ylabel("W")

plt.show()
