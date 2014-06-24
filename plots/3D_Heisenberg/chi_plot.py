#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = 4
Lx = np.arange(0, L)
Ly = np.arange(0, L)
Lz = np.arange(0, L)
r=[]
for x in Lx:
    for y in Ly:
        for z in Lz:
            r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

#Order = np.arange(1, 7)
x = 0.0

N = 64.0

Quans=["Chi"]

Files=[]

Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_1_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_2_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/3D/bare_L4_0.50_3_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

#key="Chi"
#Chi_8 = []
#for i in range(0, 6):
    #Chi_8.append(Files[i][key][0][0][0].real)

#Chi_8_pi=0.6508

#Chi_4 = []
#for i in range(6, 12):
    #Chi_4.append(Files[i][key][0][0][0].real)

#Chi_4_pi=0.6457

#for key in Quans:
    ##for i in range(len(Files)):
    ##ax.plot(Order, Files[i][key][0][0].real, label=key)
    #ax.errorbar(1.0/Order, Chi_4, marker='o', label="L=4, beta=0.90")
    #ax.errorbar(1.0/Order, Chi_8, marker='o', label="L=8, beta=0.90") 
    #ax.errorbar(0.002, Chi_4_pi, marker='*', yerr=0.0012, label="L=4, beta=0.90, Path Integral")
    #ax.errorbar(0.002, Chi_8_pi, marker='*', yerr=0.0006, label="L=8, beta=0.90, Path Integral")
#ax.set_xlim(0.0, 1.0)

#ax.errorbar(1.00, Files[12]["Chi"][0][0][0].real, marker='*', label="L=8, beta=0.90, bold")

################ Beta=0.50 ##########################
#Path=[]
#Path.append(0.716597008080860) 
#Path.append(-9.012607444665378E-002)
#Path.append(2.035853390276862E-002) 
#Path.append(-9.012607444665378E-002)

############### Beta=0.90 ##########################

#Path=[]
#Path.append(0.645719404315138)
#Path.append(-0.162600433543508)
#Path.append(7.879834033016384E-002) 
#Path.append(-0.162600433543508)

#Path=[]
#Path.append(0.65090)
#Path.append(-0.15100)
#Path.append(0.039207) 
#Path.append(-0.011242)
#Path.append(0.004776)
#Path.append(-0.011242)
#Path.append(0.039207) 
#Path.append(-0.15100)

for key in Quans:
    for i in range(len(Files)):
        ax.plot(r, Files[i][key][0].real, 'o', label=key+" Order"+str(i+1))

#ax.plot(Lx, Path, marker='*', label="Path")

ax.legend()


#plt.xlabel("1/N")
#plt.ylabel("Chi(0, 0)")

plt.xlabel("dr")
plt.ylabel("Chi(dr)")

plt.savefig("Beta0.5_L4_Chi_dr.pdf")
#plt.savefig("Chi_0_0.pdf")
plt.show()
