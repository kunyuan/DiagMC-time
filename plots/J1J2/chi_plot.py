#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = 4 
Lx = np.arange(0, L)

Order = np.arange(1, 6)

Beta = 0.50

N = 64.0

Quans=["Chi"]

Files=[]
Files.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_1_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_2_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_3_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_4_quantities.dat", Quans))
Files.append(read_data.read_array("../../data/J1J2/L4_0.50_0.50_5_quantities.dat", Quans))


fig = plt.figure()
ax = plt.subplot(111)

key = "Chi"

Chi_4 = []
for i in range(0, 5):
    Chi_4.append(Files[i][key][0][0][0].real)


#for key in Quans:
    #for i in range(len(Files)):
ax.errorbar(1.0/Order, Chi_4, marker='o', label="L=4, beta=0.50, J2=0.50")
#ax.set_xlim(0.0, 1.0)



#for key in Quans:
    #for i in range(len(Files)):
        #ax.plot(Lx, Files[i][key][0][0].real, marker='o', label=key+" Order"+str(i+1))

#ax.plot(Lx, Path, marker='*', label="Path")

ax.legend()


plt.xlabel("1/N")
plt.ylabel("Chi(0, 0)")
plt.savefig("Chi_0_0.pdf")

#plt.xlabel("dy")
#plt.ylabel("Chi(0, dy)")
#plt.savefig("Chi_dy_L4.pdf")
#plt.savefig("Chi_dy_L8.pdf")

plt.show()
