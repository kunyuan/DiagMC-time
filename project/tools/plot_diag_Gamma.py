import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)


#Gamma=["Gamma1","Gamma2","Gamma3","Gamma4"]
Gamma=["Gamma1","Gamma2","Gamma3"]
#Gamma=["Gamma1","Gamma2"]
#Gamma=["Gamma1"]

quan = read_data.read_array("./../data/bold_0.90_4_quantities.dat", Gamma)

#Gamma=["Gamma1"]
#Quans=[]
#quan = read_data.read_array("./../data/bold_0.90_1_quantities.dat")["Gamma"]
#Quans.append(quan)
#quan = read_data.read_array("./../data/bold_0.90_2_quantities.dat", Gamma)
#Quans.append(quan)
#quan = read_data.read_array("./../data/bold_0.90_3_quantities.dat", Gamma)
#Quans.append(quan)
#quan = read_data.read_array("./../data/bold_0.90_4_quantities.dat", Gamma)
#Quans.append(quan)

fig = plt.figure()
ax = plt.subplot(111)

for key in Gamma:
    ax.plot(tau, quan[key][0].diagonal().real, label=key)

#ax.plot(tau, Quans[0][0].diagonal().real, label="Order"+str(1))
#for i in range(1,4):
    #ax.plot(tau, Quans[i]["Gamma1"][0].diagonal().real, label="Order"+str(i+1))

#GamIn, _=read_data.read_array("./../0.90_Gam1.dat")["Gamma"]
#ax.plot(tau, GamIn.diagonal().real, label="NumIntGam1")

ax.legend()

plt.xlabel("tau")
plt.ylabel("diag{Gamma}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
