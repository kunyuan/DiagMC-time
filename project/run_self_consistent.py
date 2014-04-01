#!/usr/bin/python
import os
import subprocess
import time
import logging

homedir=os.getcwd()
logging.basicConfig(filename=homedir+"/project.log",level=logging.INFO,format="\n[loop.daemon][%(asctime)s][%(levelname)s]:\n%(message)s",datefmt='%y/%m/%d %H:%M:%S')
execf="./gamma3.exe"
interval=300
title = "0.50_2"
#time.sleep(300)
logging.info("Self Consistent Loop daemon started!")
logging.info(title+" is the target!")
i=0
flag = "0"

while flag=="0":
    i=i+1
    logging.info("Loop "+str(i)+" running...")
    try:
        subprocess.check_output("./collapse_data.exe")
    except subprocess.CalledProcessError as e:
        ret=e.returncode
        logging.error('collapse_data.exe return a non-zero value '+str(ret)+', something happened!')
    except:
        logging.error('Collapse_data failed!')
    else:
        try:
            os.system(execf+" <./infile/_in_selfconsist >./outfile/Output_sc")
            flag = open("stop.inp", "r").readline().lstrip().rstrip()
            #logging.info("flag: "+flag)
        except:
            logging.error('self_consistent error')

    logging.info("Loop "+str(i)+" done!")
    logging.info("Sleeping...")
    time.sleep(interval)

logging.info("Self Consistent Loop daemon ended!")


    



