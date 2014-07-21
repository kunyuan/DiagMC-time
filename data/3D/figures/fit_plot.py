import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

#Gamma=["Gamma1","Gamma2","Gamma3","Gamma4"]
#Gamma=["Gamma1","Gamma2","Gamma3"]
Gamma=["Gamma1"]
GammaBasis=["GammaBasis1"]
total=["Gamma1","GammaBasis1"]


#Gamma=["Gamma1","Gamma2"]
#Gamma=["Gamma1"]

quan = read_data.read_array("./../0.90_quantities.dat", total)

fig = plt.figure()
ax = plt.subplot(111)

for key in Gamma:
    ax.plot(tau, quan[key][0].diagonal().real, label=key)

for key in GammaBasis:
    ax.plot(tau, quan[key][0].diagonal().real, label=key)

#ax.plot(tau, Quans[0][0].diagonal().real, label="Order"+str(1))
#for i in range(1,2):
    #ax.plot(tau, Quans[i]["Gamma1"][0].diagonal().real, label=str(i))

#GamIn, _=read_data.read_array("./../0.90_Gam1.dat")["Gamma"]
#ax.plot(tau, GamIn.diagonal().real, label="NumIntGam1")

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
