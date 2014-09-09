#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 1.50
N =128

tau = np.arange(0, Beta, Beta/N)
target=["Sigma1", "Sigma2","Sigma3","Sigma4","Sigma5"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("1.50/1.50_order_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        ax.plot(tau, BoldSigma[i][key][0].real, label="{0}".format(key))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.savefig("Beta1.50_L16_J1J2_Sigma.pdf")
plt.show()