#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.80

MxT = 64
tau = np.arange(0, Beta, Beta/MxT)

fig = plt.figure()
ax = plt.subplot(111)


Quans=["ChiK"]
key = "ChiK"

########################### L=4 ########################################################
#Files=[]
#Files.append(read_data.read_array("bare_L4_0.50_1_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L4_0.50_2_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L4_0.50_3_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L4_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L4_0.50_5_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L4_0.50_6_quantities.dat", Quans))

#Order=np.arange(1, len(Files)+1)

#L = 4
#stag = (L**2+L+1)*(L/2)

#ChiL4=[]
#for i in range(0, len(Files)):
    #ChiL4.append(Files[i][key][0][stag].real)

#chi0 = 1.73477

#error = [(ChiL4[0]-chi0)*0.001, (ChiL4[1]-chi0)*0.005, (ChiL4[2]-chi0)*0.016, (ChiL4[3]-chi0)*0.012, 
        #(ChiL4[4]-chi0)*0.018, (ChiL4[5]-chi0)*0.025]
#ax.errorbar(1.0/Order, ChiL4, yerr=error, marker='o', label="L=4, beta="+str(Beta))

#######################################################################################

########################### L=8 ########################################################
Files=[]
#Files.append(read_data.read_array("bare_L8_0.50_1_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L8_0.50_2_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L8_0.50_3_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L8_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L8_0.50_5_quantities.dat", Quans))
#Files.append(read_data.read_array("bare_L8_0.50_6_quantities.dat", Quans))
Files.append(read_data.read_array("bare_L8_0.80_3_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 8
stag = (L**2+L+1)*(L/2)

ChiL8=[]
for i in range(0, len(Files)):
    ChiL8.append(Files[i][key][0][stag].real)

#chi0 = 1.74748

#error = [(ChiL8[0]-chi0)*0.001, (ChiL8[1]-chi0)*0.018, (ChiL8[2]-chi0)*0.024]
#ax.errorbar(1.0/Order, ChiL8, yerr=error, marker='^', label="L=8, beta="+str(Beta))

ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta="+str(Beta))
########################################################################################

########################### high_T ####################################################
ax.plot(0.00, 5.27, marker='*', label="beta="+str(Beta)+", high-T expansion")
#######################################################################################

ax.legend()
#ax.set_xlim(-0.05, 1.05)

plt.xlabel("1/N")
plt.ylabel("staggered susceptibility")

plt.savefig("Beta0.5_static_staggered_chi.pdf")
plt.show()
