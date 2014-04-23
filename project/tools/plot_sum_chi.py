import read_data
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Beta = 0.90
N = 64

tau = np.arange(0, Beta, Beta/N)

BareChi=[]
Chi,dim_name=read_data.read_array("./bare_0.90_2_Chi_sum.dat")["SUMChi"]
BareChi.append(Chi)
Chi,dim_name=read_data.read_array("./bare_0.90_3_Chi_sum.dat")["SUMChi"]
BareChi.append(Chi)
Chi,dim_name=read_data.read_array("./bare_0.90_4_Chi_sum.dat")["SUMChi"]
BareChi.append(Chi)

BoldChi=[]
Chi,dim_name=read_data.read_array("../data/bold_0.90_2_quantities.dat")["SUMChi"]
BoldChi.append(Chi)
Chi,dim_name=read_data.read_array("../data/bold_0.90_3_quantities.dat")["SUMChi"]
BoldChi.append(Chi)
Chi,dim_name=read_data.read_array("../data/bold_0.90_4_quantities.dat")["SUMChi"]
BoldChi.append(Chi)

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(BareChi)):
    ax.plot(tau, BareChi[i].real, label="bare"+str(i+2))

for i in range(len(BoldChi)):
    ax.plot(tau, BoldChi[i].real, label="bold"+str(i+2))

ax.legend()

plt.xlabel("tau")
plt.ylabel("sum{Chi(r)}")

# plt.savefig("0.90_Gamma.pdf")
plt.show()
