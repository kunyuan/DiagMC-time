#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

#Order = np.arange(1, 7)
x = 0.0

N = 64.0

Quans=["Chi"]
key="Chi"

fig = plt.figure()
ax = plt.subplot(111)

######################################## Chi at r=0 ####################################
#Files=[]
#Files.append(read_data.read_array("./bare_L4_0.50_1_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_2_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_3_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_5_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_6_quantities.dat", Quans))

#Chi_4 = []
#for i in range(6, 12):
    #Chi_4.append(Files[i][key][0][0][0].real)

#Chi_4_pi=0.6457

#Files=[]
#Files.append(read_data.read_array("./bare_L4_0.50_1_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_2_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_3_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_5_quantities.dat", Quans))
#Files.append(read_data.read_array("./bare_L4_0.50_6_quantities.dat", Quans))

#Chi_8 = []
#for i in range(0, 6):
    #Chi_8.append(Files[i][key][0][0][0].real)

#Chi_8_pi=0.6508



#for key in Quans:
    #ax.errorbar(1.0/Order, Chi_4, marker='o', label="L=4, beta=0.90")
    #ax.errorbar(1.0/Order, Chi_8, marker='o', label="L=8, beta=0.90") 
    #ax.errorbar(0.0, Chi_4_pi, marker='*', yerr=0.0012, label="L=4, beta=0.90, Path Integral")
    #ax.errorbar(0.0, Chi_8_pi, marker='*', yerr=0.0006, label="L=8, beta=0.90, Path Integral")
#ax.set_xlim(-0.05, 1.05)
#ax.legend()

#plt.xlabel("1/N")
#plt.ylabel("Chi(0, 0)")

#plt.savefig("Chi_0_0.pdf")
################################################################################################

########################################L=16 Chi(r)############################################
Files=[]
Files.append(read_data.read_array("./L16_0.70_1_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_2_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_3_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_4_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L16_0.70_5_bare_quantities.dat", Quans))

L = 16
Lx = np.arange(0, L)
Ly = np.arange(0, L)
r=[]
for x in Lx:
    for y in Ly:
        r.append(np.sqrt(x**2.0+y**2.0))

for key in Quans:
    for i in range(len(Files)):
        ax.plot(r, Files[i][key][0].real, 'o', label=key+" Order"+str(i+1))
ax.legend()

plt.xlabel("dr")
plt.ylabel("Chi(dr)")

plt.savefig("Beta0.7_L16_Chi_dr.pdf")
######################################################################################



plt.show()

