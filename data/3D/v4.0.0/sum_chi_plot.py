#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.80
N = 64

tau = np.arange(0, Beta, Beta/N)

fig = plt.figure()
ax = plt.subplot(111)

Quans=["SUMChi"]
key = "SUMChi"

SUMChi=[]
SUMChi.append(read_data.read_array("./bare_L8_0.80_3_quantities.dat", Quans))

for i in range(len(SUMChi)):
    ax.plot(tau, SUMChi[i][key][0].real, marker='o', label="L=8, beta="+str(Beta))

ax.legend()

plt.xlabel("1/N")
plt.ylabel("Chi")

plt.savefig("Beta"+str(Beta)+"_static_uniform_chi.pdf")
plt.show()
