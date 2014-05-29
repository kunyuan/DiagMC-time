#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

L = np.arange(0, 4)
N = 64

Quans=["Chi"]

Files=[]
#Files.append(read_data.read_array("./0.90_2_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_4_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_5_bold_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.90_4_bare_quantities.dat", Quans))

#Files.append(read_data.read_array("./0.50_0.10_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.50_0.50_4_quantities.dat", Quans))
#Files.append(read_data.read_array("./0.50_1.00_4_quantities.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(Files)):
    for key in Quans:
        #ax.plot(L, Files[i][key][0][0].real, label=key)
        chi=[]
        for j in range(0, 4):
            chi.append(sum(Files[i][key][0][0][j]).real/N)
        ax.plot(L, chi, label=key)

################ Beta=0.50 ##########################
#Path=[]
#Path.append(0.716597008080860) 
#Path.append(-9.012607444665378E-002)
#Path.append(2.035853390276862E-002) 
#Path.append(-9.012607444665378E-002)

################ Beta=0.90 ##########################
Path=[]
Path.append(0.643467765710929 )
Path.append(-0.164205140020811)
Path.append(8.183978106942874E-002) 
Path.append(-0.164205140020811)

ax.plot(L, Path, label="Path")

ax.legend()

plt.xlabel("dy")
plt.ylabel("Chi(0, dy)")

#plt.savefig("0.90_bold_Chi.pdf")
#plt.savefig("0.90_bare_Chi.pdf")
plt.show()
