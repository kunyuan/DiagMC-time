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

Quans=["Polar"]

Polar=read_data.read_array("../0.90_quantities.dat", Quans)

fig = plt.figure()
ax = plt.subplot(111)

for key in Quans:
    ax.plot(tau, Polar[key][0][0,0,:].real, label=key)
    ax.plot(tau, Polar[key][0][0,0,:].imag, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("Pi(0,0)")

plt.show()
