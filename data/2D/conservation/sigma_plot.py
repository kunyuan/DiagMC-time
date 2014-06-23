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

BoldSigma=[]
Sigma,dim_name=read_data.read_array("./0.90_2_quantities.dat")["Sigma"]
BoldSigma.append(Sigma)
Sigma,dim_name=read_data.read_array("./0.90_3_quantities.dat")["Sigma"]
BoldSigma.append(Sigma)
Sigma,dim_name=read_data.read_array("./0.90_4_quantities.dat")["Sigma"]
BoldSigma.append(Sigma)

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldSigma)):
    ax.plot(tau, BoldSigma[i].real, label="Order"+str(i+2))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
