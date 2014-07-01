#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

def get_r(x,y,z,L):
    r = z*L**2+y*L+x
    return r


x = 0.0

N = 64.0

Quans=["GammaR"]
key = "GammaR"


##################################L=4#################################
#L = 4
#Lx = np.arange(0, L)
#Ly = np.arange(0, L)
#Lz = np.arange(0, L)

#r=[]
#for x in Lx:
    #for y in Ly:
        #for z in Lz:
            #r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

#Files=[]
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.70_1_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.70_2_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.70_3_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.70_4_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.70_5_quantities.dat", Quans))
#Files.append(read_data.read_array("../../data/3D/bare_L4_0.70_6_quantities.dat", Quans))

#fig = plt.figure() 
#ax = plt.subplot(111)

#for key in Quans:
    #for i in range(len(Files)):
        #ax.plot(r, Files[i][key][0].real, 'o', label=key+" Order"+str(i+1))
#ax.legend()

#plt.xlabel("dr") 
#plt.ylabel("Gamma(dr)")

#plt.savefig("Beta0.7_L4_Gamma_r.pdf")
########################################################################


##################################L=8#################################
L = 8
Lx = np.arange(0, L/2+1)
Ly = np.arange(0, L/2+1)
Lz = np.arange(0, L/2+1)

Files=[]
Files.append(read_data.read_array("../../../data/3D/bare_L8_0.70_1_quantities.dat", Quans))
Files.append(read_data.read_array("../../../data/3D/bare_L8_0.70_2_quantities.dat", Quans))
Files.append(read_data.read_array("../../../data/3D/bare_L8_0.70_3_quantities.dat", Quans))
Files.append(read_data.read_array("../../../data/3D/bare_L8_0.70_4_quantities.dat", Quans))
Files.append(read_data.read_array("../../../data/3D/bare_L8_0.70_5_quantities.dat", Quans))

Gamma=[]
for i in range(len(Files)):
    r=[]
    Gammai=[]
    for x in Lx:
        for y in np.arange(0,x+1):
            for z in np.arange(0,y+1):
                r.append(np.sqrt(x**2.0+y**2.0+z**2.0))
                Gammai.append(Files[i][key][0][get_r(x,y,z,L)].real)
    Gamma.append(Gammai)

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(Files)):
    ax.plot(r, Gamma[i], 'o', label=key+" Order"+str(i+1))

ax.set_yscale('log')
ax.legend()


plt.xlabel("dr") 
plt.ylabel("Gamma(dr)")

plt.savefig("Beta0.7_L8_Gamma_r.pdf")
########################################################################



plt.show()

