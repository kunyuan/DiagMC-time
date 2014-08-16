#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

def get_r(x,y,z,L):
    r = z*L**2+y*L+x
    return r

N = 64.0

Quans=["Chi"]
key="Chi"

fig = plt.figure()
ax = plt.subplot(111)

Files=[]
Files.append(read_data.read_array("bare_L8_0.80_3_quantities.dat", Quans))

L = 8
Lx = np.arange(0, L/2)
Ly = np.arange(0, L/2)
Lz = np.arange(0, L/2)
r=[]
realChi =[]
imagChi =[]
for x in Lx:
    for y in Ly:
        for z in Lz:
            r.append(np.sqrt(x**2.0+y**2.0+z**2.0))
            realChi.append(Files[0][key][0][get_r(x,y,z,L)].real)
            imagChi.append(Files[0][key][0][get_r(x,y,z,L)].imag)

ax.plot(r, realChi, 'o')
ax.legend()

plt.xlabel("dr")
plt.ylabel("Chi(dr)")

plt.savefig("Beta0.8_L8_Chi_dr.pdf")

plt.show()

################ Beta=0.50 ##########################
#Path=[]
#Path.append(0.716597008080860) 
#Path.append(-9.012607444665378E-002)
#Path.append(2.035853390276862E-002) 
#Path.append(-9.012607444665378E-002)

############### Beta=0.90 ##########################

#Path=[]
#Path.append(0.645719404315138)
#Path.append(-0.162600433543508)
#Path.append(7.879834033016384E-002) 
#Path.append(-0.162600433543508)

#Path=[]
#Path.append(0.65090)
#Path.append(-0.15100)
#Path.append(0.039207) 
#Path.append(-0.011242)
#Path.append(0.004776)
#Path.append(-0.011242)
#Path.append(0.039207) 
#Path.append(-0.15100)
