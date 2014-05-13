#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

Gamma=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Gamma=["Gamma1","Gamma2","Gamma3"]
#Gamma=["Gamma1","Gamma2","Gamma"]
#Gamma=["Gamma"]

#quan = read_data.read_array("./../../data/conservation/bare_0.90_4_quantities.dat", Gamma)
quan = read_data.read_array("./bare_2/0.90_quantities.dat", Gamma)
#quan = read_data.read_array("../0.90_Gam1.dat", Gamma)
#quan = read_data.read_array("../0.90_quantities.dat", Gamma)

#Gamma=["Gamma1"]
#Quans=[]
#quan = read_data.read_array("./../data/bold_0.90_1_quantities.dat")["Gamma"]
#Quans.append(quan)
#quan = read_data.read_array("./../data/bold_0.90_2_quantities.dat", Gamma)
#Quans.append(quan)
#quan = read_data.read_array("./../data/bold_0.90_3_quantities.dat", Gamma)
#Quans.append(quan)
#quan = read_data.read_array("./../0.90_quantities.dat", Gamma)
#Quans.append(quan)
#quan = read_data.read_array("./../data/bold_0.90_4_quantities.dat", Gamma)
#Quans.append(quan)

fig = plt.figure()
ax = plt.subplot(111)

for key in Gamma:
    #ax.plot(tau, quan[key][0].diagonal().imag, label=key)
    ax.plot(tau, quan[key][0].diagonal().real, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
