#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

#Order = 1
#with open("bare_L8_0.80_1_denom.dat") as f:
    #data = f.read()

#Order = 2
#with open("bare_L8_0.80_2_denom.dat") as f:
    #data = f.read()

Order = 3
with open("bare_L8_0.80_3_denom.dat") as f:
    data = f.read()
data = data.split('\n')
real=[]
for i in range(0, len(data)-1):
    real.append(float(data[i].split('      ')[0]))


#with open("bare_L8_0.80_3_0.01_denom.dat") as f:
    #data2 = f.read()
#data2 = data2.split('\n')
#real2=[]
#for i in range(0, len(data2)-1):
    #real2.append(float(data2[i].split('      ')[0]))

fig = plt.figure()
ax = plt.subplot(111)


ax.plot(range(0, len(real)), real, marker='o', label="L=8, Beta=(0.65, 0.05, 0.80), Order="+str(Order))
#ax.plot(range(0, len(real2)), real2, marker='o', label="L=8, Beta=(0.65, 0.01, 0.80), Order="+str(Order))

ax.legend()

plt.xlabel("time")
plt.ylabel("Denom(pi, pi, pi)")

plt.savefig("Beta0.8_L8_"+str(Order)+"_denom_time.pdf")
plt.show()
