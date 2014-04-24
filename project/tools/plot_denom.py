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

BoldDenom=[]
Denom,dim_name=read_data.read_array("../data/bold_0.90_2_quantities.dat")["Denom"]
BoldDenom.append(Denom[:][L/2][L/2])
Denom,dim_name=read_data.read_array("../data/bold_0.90_3_quantities.dat")["Denom"]
BoldDenom.append(Denom[:][L/2][L/2])
Denom,dim_name=read_data.read_array("../data/bold_0.90_4_quantities.dat")["Denom"]
BoldDenom.append(Denom[:][L/2][L/2])

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BoldDenom)):
    ax.plot(tau, BoldDenom[i].real, label="Order"+str(i+2))

ax.legend()

plt.xlabel("tau")
plt.ylabel("Sigma")

plt.show()
