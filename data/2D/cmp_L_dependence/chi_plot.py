#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = 8
Lx = np.arange(0, L)

Order = np.arange(1, 6)
x = 0.0

N = 64.0

Quans=["Chi"]

Files=[]
#Files.append(read_data.read_array("./0.90_2_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_4_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_5_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_4_bare_quantities.dat", Quans))

#Files.append(read_data.read_array("./0.50_0.10_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.50_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("../L8_0.90_4_bare_quantities.dat", Quans))
Files.append(read_data.read_array("./L8_0.90_1_quantities.dat", Quans))
Files.append(read_data.read_array("./L8_0.90_2_quantities.dat", Quans))
Files.append(read_data.read_array("./L8_0.90_3_quantities.dat", Quans))
Files.append(read_data.read_array("./L8_0.90_4_quantities.dat", Quans))
Files.append(read_data.read_array("./L8_0.90_5_quantities.dat", Quans))

#Files.append(read_data.read_array("./L4_0.90_1_quantities.dat", Quans))
#Files.append(read_data.read_array("./L4_0.90_2_quantities.dat", Quans))
#Files.append(read_data.read_array("./L4_0.90_3_quantities.dat", Quans))
#Files.append(read_data.read_array("./L4_0.90_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./L4_0.90_5_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

#key="Chi"
#Chi_8 = []
#for i in range(0, 5):
    #Chi_8.append(Files[i][key][0][0][0].real)

#Chi_8_pi=0.6508

#Chi_4 = []
#for i in range(5, 10):
    #Chi_4.append(Files[i][key][0][0][0].real)

#Chi_4_pi=0.6435

#for key in Quans:
    ##for i in range(len(Files)):
    ##ax.plot(Order, Files[i][key][0][0].real, label=key)
    #ax.errorbar(1.0/Order, Chi_4, marker='o', label="L=4, beta=0.90")
    #ax.errorbar(1.0/Order, Chi_8, marker='o', label="L=8, beta=0.90") 
    #ax.errorbar(0.2, Chi_4_pi, marker='*', yerr=0.0012, label="L=4, beta=0.90, Path Integral")
    #ax.errorbar(0.2, Chi_8_pi, marker='*', yerr=0.0006, label="L=8, beta=0.90, Path Integral")
#ax.set_xlim(0.0, 1.0)

################ Beta=0.50 ##########################
#Path=[]
#Path.append(0.716597008080860) 
#Path.append(-9.012607444665378E-002)
#Path.append(2.035853390276862E-002) 
#Path.append(-9.012607444665378E-002)

############### Beta=0.90 ##########################

#Path=[]
#Path.append(0.643467765710929)
#Path.append(-0.164205140020811)
#Path.append(8.183978106942874E-002) 
#Path.append(-0.164205140020811)

Path=[]
Path.append(0.65090)
Path.append(-0.15100)
Path.append(0.039207) 
Path.append(-0.011242)
Path.append(0.004776)
Path.append(-0.011242)
Path.append(0.039207) 
Path.append(-0.15100)

for key in Quans:
    for i in range(len(Files)):
        ax.plot(Lx, Files[i][key][0][0].real, marker='o', label=key+" Order"+str(i+1))

ax.plot(Lx, Path, marker='*', label="Path")


ax.legend()

#plt.xlabel("1/N")
#plt.ylabel("Chi(0, 0)")

plt.xlabel("dy")
plt.ylabel("Chi(0, dy)")

#plt.savefig("Chi_dy_L4.pdf")
plt.savefig("Chi_dy_L8.pdf")
#plt.savefig("Chi_0_0.pdf")
#plt.show()
