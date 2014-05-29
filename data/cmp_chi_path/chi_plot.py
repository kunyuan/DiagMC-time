#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import read_data

Beta = 0.90
N = 64

L = np.arange(0, 4)

Quans=["Chi1", "Chi2", "Chi3", "Chi4"]

Files=[]
Files.append(read_data.read_array("./bold_Chi.dat", Quans))
#Files.append(read_data.read_array("./bare_Chi.dat", Quans))

fig = plt.figure()
ax = plt.subplot(111)

for i in range(len(Files)):
    for key in Quans:
        ax.plot(L, Files[i][key][0][0].real/Beta, label=key)

Path=[]
Path.append(0.643467765710929) 
Path.append(-0.164205140020811)      
Path.append(8.183978106942874e-002) 
Path.append(-0.164205140020811)      

ax.plot(L, Path, label="Path-Integral")

ax.legend()

plt.xlabel("dy")
plt.ylabel("Chi(0, dy)")

#plt.savefig("0.90_bold_Chi.pdf")
plt.savefig("0.90_bare_Chi.pdf")
plt.show()
