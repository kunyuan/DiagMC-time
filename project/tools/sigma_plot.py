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
target=["Sigma"]

BoldSigma=[]
BoldSigma.append(read_data.read_array("./../0.90_quantities.dat",target))
BoldSigma.append(read_data.read_array("./../old_0.90_quantities.dat",target))
#BoldSigma.append(read_data.read_array("./../bold_0.90_1_quantities.dat",target))
#BoldSigma.append(read_data.read_array("./../bold_0.90_2_quantities.dat",target))
#BoldSigma.append(read_data.read_array("./../bold_0.90_3_quantities.dat",target))
#BoldSigma.append(read_data.read_array("./0.90_4_bare_quantities.dat",target))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    for key in target:
        ax.plot(tau, BoldSigma[i][key][0].real, label="{0}, Order{1}".format(key, i))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
