import os
#title = raw_input("Enter the file to collapse:")
title = "0.50_2"
print title+" is the target!"
file1 = title + "_Gamma_MC.dat"
file2 = title + "_-*_Gamma_MC.dat"
file3 = title + "_monte_carlo_data.dat"
file4 = title + "_-*_monte_carlo_data.dat"

os.system("rm data_collapse/"+file1+"\n")
os.system("cat "+file2+" >> data_collapse/"+file1+"\n")
os.system("rm data_collapse/"+file3+"\n")
os.system("cat "+file4+" >> data_collapse/"+file3+"\n")
os.system("echo "+title+"|./collapse_data")


    



