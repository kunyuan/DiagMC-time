#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.90
N = 64
Quans=["Sigma"]

tau = np.arange(0, Beta, Beta/N)
tau2 = np.arange(Beta*1.0/N, Beta+Beta*1.0/N, Beta/N)

BoldSigma=[]
BoldSigma.append(read_data.read_array("../0.90_quantities.dat",Quans))
#BoldSigma.append(read_data.read_array("./0.90_before.dat",Quans))

fig = plt.figure()
ax = plt.subplot(111)

#key="Sigma"
#ax.plot(tau, BoldSigma[key][0].real, label=key)
#key="G"
#ax.plot(tau, BoldSigma[key][0].real*0.531527, label=key)

for i in range(len(BoldSigma)):
    for key in Quans:
        ax.plot(tau, BoldSigma[i][key][0].imag, label=key+str(i))
        ax.plot(tau, BoldSigma[i][key][0].real, label=key+str(i))
        #ax.plot(tau2, BoldSigma[i][key][0][::-1].real, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
