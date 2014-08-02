#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.50

MxT = 64
tau = np.arange(0, Beta, Beta/MxT)

Quans=["ChiKt0"]
key = "ChiKt0"

Files=[]
Files.append(read_data.read_array("bare_L8_0.50_3_quantities.dat", Quans))
Files.append(read_data.read_array("bare_L8_0.50_4_quantities.dat", Quans))
Files.append(read_data.read_array("bare_L8_0.50_5_quantities.dat", Quans))

Order=np.arange(1, len(Files)+1)

L = 8
stag = (L**2+L+1)*(L/2)

ChiL8=[]
for i in range(0, len(Files)):
    ChiL8.append(Files[i][key][0][stag].real)
    print i, ChiL8[i]*L**3.0/3.0

#fig = plt.figure() 
#ax = plt.subplot(111)


#ax.plot(1.0/Order, ChiL8, marker='o', label="L=8, beta="+str(Beta))
########################################################################################

#ax.legend()

#plt.xlabel("1/N")
#plt.ylabel("staggered susceptibility")

#plt.savefig("Beta0.5_static_staggered_chi.pdf")
#plt.show()
