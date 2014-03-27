import os
title = raw_input("Enter the file to collapse:")
#title = "0.50_2"

print title+" is the target!"

file3 = title + "_monte_carlo_data.bin.dat"
file4 = "*_" + title + "_monte_carlo_data.bin.dat"

os.system("rm "+file3+"\n")
os.system("cat ../"+file4+" >> "+file3+"\n")
os.system("./collapse_data.exe")


    



