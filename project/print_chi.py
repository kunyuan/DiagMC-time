#!/usr/bin/env python
import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

L = 8
Beta = 0.70
MxT = 128
tau = np.arange(0, Beta, Beta/MxT)

MnOrder = 2
MxOrder = 5
Quans=["ChiK"]
key="ChiK"

stag = (L**2+L+1)*(L/2)
unif = 0

Files=[]
for i in range(MnOrder, MxOrder+1):
    Files.append(read_data.read_array("L"+str(L)+"_"+"{:.2f}".format(Beta)+"_"+str(i)+"/"+"{:.2f}".format(Beta)+"_quantities.dat", Quans))
    print Files[i][key][0][stag]
    print Files[i][key][0][unif]





########################################################################################

