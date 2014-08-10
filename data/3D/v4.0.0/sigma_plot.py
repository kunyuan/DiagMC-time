#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.70
N =128

tau1 = np.arange(0, Beta/2, Beta/N)
tau2 = np.arange(Beta/2, Beta, Beta/N)
target=["Sigma"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("bare_L8_0.70_5_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        #ax.plot(tau, BoldSigma[i][key][0].real, label="{0}, Order{1}".format(key, i+1))
        ax.plot(tau1, -1.0*BoldSigma[i][key][0][0:N/2].imag, label="{0}, Order{1}".format(key, i+1))
        ax.plot(tau2, BoldSigma[i][key][0][N/2:N].imag, label="{0}, Order{1}".format(key, i+1))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.savefig("Beta0.70_L8_Sigma.pdf")
plt.show()
