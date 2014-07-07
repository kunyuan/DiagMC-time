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
target=["Sigma"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("../../../../project/0.50_quantities.dat",target))
BoldSigma.append(read_data.read_array("../../../../old_program/0.50_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        ax.plot(tau, BoldSigma[i][key][0].real, label="{0}, Order{1}".format(key, i+1))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.savefig("Beta0.5_L4_Sigma.pdf")
plt.show()
