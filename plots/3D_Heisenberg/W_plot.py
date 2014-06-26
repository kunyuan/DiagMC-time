#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.70
N = 64

tau = np.arange(0, Beta, Beta/N)
target=["W"]

fig = plt.figure()
ax = plt.subplot(111)


############################# L=4 ##################################################
BoldSigma=[]
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_1_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_2_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_3_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_4_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_5_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.50_6_quantities.dat",target))
BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.70_1_quantities.dat",target))
BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.70_2_quantities.dat",target))
BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.70_3_quantities.dat",target))
BoldSigma.append(read_data.read_array("../../data/3D/bare_L4_0.70_4_quantities.dat",target))

L = 4
Lx = range(0, L)
Ly = range(0, L)
Lz = range(0, L)
Vol = L*L*L

W=[]
key = "W"
for i in range(len(BoldSigma)):
    Wr=[]
    for r in range(0, Vol):
        Wr.append(abs(sum(BoldSigma[i][key][0][r]).real))
    W.append(Wr)

r=[]
for x in Lx:
    for y in Ly:
        for z in Lz:
            r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

for i in range(len(BoldSigma)):
    ax.plot(r, W[i], 'o', label="{0}, Beta={1}, Order{2}".format(key, Beta, i+1))

ax.legend()

plt.xlabel("dr")
plt.ylabel("|W(dr)|")
plt.savefig("Beta"+str(Beta)+"_L4_W_r.pdf")
################################################################################################



############################# L=8 ##################################################
#BoldSigma=[]
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L8_0.50_1_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L8_0.50_2_quantities.dat",target))
#BoldSigma.append(read_data.read_array("../../data/3D/bare_L8_0.50_3_quantities.dat",target))

#L = 8
#Lx = range(0, L)
#Ly = range(0, L)
#Lz = range(0, L)
#Vol = L*L*L

#W=[]
#key = "W"
#for i in range(len(BoldSigma)):
    #Wr=[]
    #for r in range(0, Vol):
        #Wr.append(abs(sum(BoldSigma[i][key][0][r]).real))
    #W.append(Wr)

#r=[]
#for x in Lx:
    #for y in Ly:
        #for z in Lz:
            #r.append(np.sqrt(x**2.0+y**2.0+z**2.0))

#for i in range(len(BoldSigma)):
    #ax.plot(r, W[i], 'o', label="{0}, Beta=0.5, Order{1}".format(key, i+1))

#ax.legend()

#plt.xlabel("dr")
#plt.ylabel("|W(dr)|")
#plt.savefig("Beta0.5_L8_W_r.pdf")
################################################################################################
plt.show()
