#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = plt.subplot(111)

for Order in range(3, 6):
    data = np.loadtxt("bare_L8_0.50_"+str(Order)+"_denom.dat")

    real=[]
    for i in range(0, len(data)-1):
        real.append(data[i][0])

    ax.plot(range(0, len(real)), real, marker='o', label="L=8, Beta=0.50, Order="+str(Order))

ax.legend()

plt.xlabel("time")
plt.ylabel("Denom(pi, pi, pi)")

plt.savefig("Beta0.5_L8_denom_time.pdf")
plt.show()
