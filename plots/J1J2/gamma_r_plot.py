#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

x = 0.0

N = 64.0

Quans=["GammaR"]


fig = plt.figure()
ax = plt.subplot(111)

##################################L=4#################################
#L = 4
#Lx = np.arange(0, L)
#Ly = np.arange(0, L)
#Lz = np.arange(0, L)

#r=[]
#for x in Lx:
    #for y in Ly:
        #for z in Lz:
            #r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

#Files=[]
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_1_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_2_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_3_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_5_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_6_quantities.dat", Quans))
##Files.append(read_data.read_array("../../project/0.70_quantities.dat", Quans))

#for key in Quans:
    #for i in range(len(Files)):
        #ax.plot(r, Files[i][key][0].real, 'o', label=key+" Order"+str(i+1))
#ax.legend()

#plt.xlabel("dr") 
#plt.ylabel("Gamma(dr)")

#plt.savefig("Beta0.5_L4_Gamma_r.pdf")
########################################################################


##################################L=8#################################
L = 8
Lx = np.arange(0, L)
Ly = np.arange(0, L)
Lz = np.arange(0, L)

r=[]
for x in Lx:
    for y in Ly:
        for z in Lz:
            r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

Files=[]
Files.append(read_data.read_array("../../data/3D/bare_L8_0.50_1_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L8_0.50_2_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L8_0.50_3_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)


for key in Quans:
    for i in range(len(Files)):
        ax.plot(r, Files[i][key][0].real, 'o', label=key+" Order"+str(i+1))
ax.legend()


plt.xlabel("dr") 
plt.ylabel("Gamma(dr)")

plt.savefig("Beta0.5_L8_Gamma_r.pdf")
########################################################################



plt.show()
