import os
import time
import logging

logging.basicConfig(filename = os.path.join(os.getcwd(), 'self_consistent.log'),filemode = 'w', format = '%(asctime)s - %(levelname)s: %(message)s')

execf="./gamma3.exe"
interval=300
title = "0.50_2"
#time.sleep(300)
logging.info(title+" is the target!")
i=0
flag = "0"

while flag=="0":
    try:
        os.system("./collapse_data.exe")
    except:
        logging.error('collapse_data error')

    try:
        os.system(execf+" <./infile/_in_selfconsist >./outfile/Output_sc")
        flag = open("stop.inp", "r").readline().lstrip().rstrip()
        logging.info("flag: "+flag)
    except:
        logging.error('self_consistent error')
    time.sleep(interval)

logging.info("self_consistent loop is stopped!")


    



