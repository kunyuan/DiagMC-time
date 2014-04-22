#!/usr/bin/python
'''This code could be used to run self consistent loop or output loop'''
import os
import sys
import subprocess
import time
import logging

INTERVAL = 300
EXEC = "./gamma3.exe"

def run_loop(infile, outfile):
    '''the loop to do self consisent calculation or output variables'''
    homedir = os.getcwd()
    logging.basicConfig(filename=homedir+"/project.log", level=logging.INFO,
         format="\n[loop.daemon][%(asctime)s][%(levelname)s]:\n%(message)s",
         datefmt='%y/%m/%d %H:%M:%S')
    execf = os.path.abspath(EXEC)
    logging.info("Loop daemon started!")
    #logging.info(title+" is the target!")
    print "loop daemon started..."
    i = 0
    while True:
        i += 1
        logging.info("Loop "+str(i)+" running...")
        try:
            proc = subprocess.Popen(homedir+"/collapse_data.exe")
            exitcode = proc.wait()
            if exitcode is not 0:
                raise RuntimeError('collapse_data return a value '+exitcode)
        except subprocess.CalledProcessError as err:
            ret = err.returncode
            logging.error('collapse_data.exe return a non-zero value '
                +str(ret)+', something happened!')
        except RuntimeError as err:
            logging.error(err.message)
        except:
            logging.error('Collapse_data failed!')
        else:
            try:
                #os.system('echo '+infile+' | '+execf+" >"+outfile)
                os.system("exec "+execf+"<"+infile+" >>"+outfile)
            except:
                logging.error('loop error!')

        logging.info("Loop "+str(i)+" done!")
        logging.info("Sleeping...")
        time.sleep(INTERVAL)
        logging.info("Loop daemon ended!")

INFILE = raw_input("Please input the infile path: ")
INFILE = os.path.abspath(INFILE)
TYPENAME = INFILE.split("/")[-1]
if TYPENAME.find('_in') is -1:
    if TYPENAME.find('in') is -1: 
        OUTFILE = "&1"
    else:
        TYPENAME = TYPENAME.replace('in','out')+'.txt'
        OUTFILE = '/'.join(INFILE.split('/')[:-2])+'/'+TYPENAME
else:
    TYPENAME = TYPENAME.replace('_in','out')+'.txt'
    OUTFILE = '/'.join(INFILE.split('/')[:-2])+'/'+TYPENAME

print INFILE, OUTFILE

print "Start loop..."
run_loop(INFILE, OUTFILE)
