#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

#Gamma=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Gamma=["Gamma1","Gamma2","Gamma3"]
#Gamma=["Gamma1","Gamma2"]
#Gamma=["Gamma", "Gamma1", "GammaBasis1"]
Gamma=["Gamma", "Gamma1"]
#Gamma=["Gam0PF0","Gam0PF1","Gam0PF2"]

#quan = read_data.read_array("./../../data/conservation/bare_0.90_4_quantities.dat", Gamma)
#quan = read_data.read_array("../bold_0.90_2_quantities.dat", Gamma)
#quan = read_data.read_array("../bold_0.90_4_quantities.dat", Gamma)
#quan = read_data.read_array("./bold_0.90_4_quan.dat", Gamma)
quan = read_data.read_array("../0.90_quantities.dat", Gamma)
#quan = read_data.read_array("./0.90_quantities.dat", Gamma)
#quan = read_data.read_array("../program/test/test_fourier.dat", Gamma)


fig = plt.figure()
ax = plt.subplot(111)

for key in Gamma:
    ax.plot(tau, quan[key].data().diagonal().real, label=key)

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
