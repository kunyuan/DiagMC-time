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

Quans=["SUMChi"]

BoldChi=[]
BoldChi.append(read_data.read_array("./0.90_2_bold_quantities.dat", Quans))
BoldChi.append(read_data.read_array("./0.90_4_bold_quantities.dat", Quans))
BoldChi.append(read_data.read_array("./0.90_4_bare_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldChi)):
    for key in Quans:
        ax.plot(tau, BoldChi[i][key][0].real, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("sum{Chi(r)}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
