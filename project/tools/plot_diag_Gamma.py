import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

#Gamma=["Gamma1","Gamma2"]
Gamma=["Gamma1"]
quan = read_data.read_array("./../0.90_quantities.dat", Gamma)

fig = plt.figure()
ax = plt.subplot(111)

for key in Gamma:
    ax.plot(tau, quan[key][0].diagonal().real, label=key)
GamIn, _=read_data.read_array("./../0.90_Gam1.dat")["Gamma"]
ax.plot(tau, GamIn.diagonal().real, label="NumIntGam1")

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
