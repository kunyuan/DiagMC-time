import os
import time
execf="./gamma3"
interval=1200
title = "0.50_2"
#time.sleep(3600)
print title+" is the target!"
file1 = title + "_Gamma_MC.dat"
file2 = title + "_-*_Gamma_MC.dat"
file3 = title + "_monte_carlo_data.dat"
file4 = title + "_-*_monte_carlo_data.dat"
i=0
flag = "0"
while flag=="0":
    os.system("rm data_collapse/"+file1+"\n")
    os.system("cat "+file2+" >> data_collapse/"+file1+"\n")
    os.system("rm data_collapse/"+file3+"\n")
    os.system("cat "+file4+" >> data_collapse/"+file3+"\n")
    os.system("echo "+title+"|./collapse_data")
    os.system(execf+" <./infile/_in_selfconsist >Output_sc")
    flag = open("stop", "r").readline().lstrip().rstrip()
    print flag
    time.sleep(interval)
print "Self consistent loop is stopped!"


    



