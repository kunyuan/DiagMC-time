import os
import time
execf="./gamma3.exe"
interval=300
title = "0.50_2"
#time.sleep(3600)
print title+" is the target!"
i=0
flag = "0"

file3 = title + "_monte_carlo_data.bin.dat"
file4 = "*_" + title + "_monte_carlo_data.bin.dat"

while flag=="0":
    os.system("rm data_process/"+file3+"\n")
    os.system("cat "+file4+" >> data_process/"+file3+"\n")
    os.system("echo "+title+"|./data_process/collapse_data.exe")

    os.system(execf+" <./infile/_in_selfconsist >Output_sc")
    flag = open("stop", "r").readline().lstrip().rstrip()
    print flag
    time.sleep(interval)
print "Self consistent loop is stopped!"


    



