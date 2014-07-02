#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.70

MxT = 64
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)


Quans=["ChiK"]

########################### L=4 ########################################################
Files=[]

Files.append(read_data.read_array("./bare_L4_0.70_1_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L4_0.70_2_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L4_0.70_3_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L4_0.70_4_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L4_0.70_5_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L4_0.70_6_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 4
stag = (L**2+L+1)*(L/2)

key = "ChiK"
ChiL4=[]
for i in range(0, len(Files)):
    ChiL4.append(Files[i][key][0][stag].real)

ax.plot(1.0/Order, ChiL4, marker='o', label="L=4, beta="+str(Beta))
#######################################################################################

########################### L=8 ########################################################
Files=[]
Files.append(read_data.read_array("./bare_L8_0.70_1_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L8_0.70_2_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L8_0.70_3_quantities.dat", Quans))
Files.append(read_data.read_array("./bare_L8_0.70_4_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 8
stag = (L**2+L+1)*(L/2)

ChiL8=[]
for i in range(0, len(Files)):
    ChiL8.append(Files[i][key][0][stag].real)

ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta="+str(Beta))
########################################################################################

########################### high_T ####################################################
ax.plot(0.00, 3.65, marker='*', label="beta="+str(Beta)+", high-T expansion")
#######################################################################################

ax.legend()
ax.set_xlim(-0.05, 1.05)

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.7_static_staggered_chi.pdf")
plt.show()
