#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.70
N =128

tau = np.arange(Beta/(2*N), Beta+Beta/(2*N), Beta/N)
#target=["Sigma1", "Sigma2","Sigma3","Sigma4","Sigma5","Sigma6"]
target=["Sigma2","Sigma4"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("L8_0.70_6/0.70_order_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        ax.plot(tau, BoldSigma[i][key][0].real, label="{0}".format(key))
        ax.plot(tau, BoldSigma[i][key][0].imag, label="{0}".format(key))
        ax.plot(tau, BoldSigma[i][key][0][::-1].real, label="{0}".format(key))
        ax.plot(tau, -BoldSigma[i][key][0][::-1].imag, label="{0}".format(key))
        #ax.plot(tau, (BoldSigma[i][key][0].real-BoldSigma[i][key][0][::-1].real)/BoldSigma[i][key][0][N/2].real, label="{0}".format(key))
        #ax.plot(tau, (BoldSigma[i][key][0].imag+BoldSigma[i][key][0][::-1].imag)/BoldSigma[i][key][0][0].imag, label="{0}".format(key))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.savefig("Beta0.70_L8_Order24_Sigma.pdf")
plt.show()
