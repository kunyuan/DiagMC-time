#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.90
N = 64
L = 4

tau = np.arange(0, Beta, Beta/N)
target=["W"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("../0.90_quantities.dat",target))
#BoldSigma.append(read_data.read_array("./0.90_2_quantities.dat",target))
#BoldSigma.append(read_data.read_array("./0.90_3_quantities.dat",target))
#BoldSigma.append(read_data.read_array("./0.90_4_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        print BoldSigma[i][key][0].shape
        ax.plot(tau, BoldSigma[i][key][0][0,0,:].real, label="{0}, Order{1}".format(key, i+2))
        ax.plot(tau, BoldSigma[i][key][0][1,0,:].real, label="{0}, Order{1}".format(key, i+2))
        ax.plot(tau, BoldSigma[i][key][0][1,1,:].real, label="{0}, Order{1}".format(key, i+2))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
