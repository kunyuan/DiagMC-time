#!/usr/bin/python
'''This is the job manage script.
   This is a quite universal code for all different type of simulations'''
import os
import time
import subprocess
import logging

import inlist

PROCLIST = []
logging.basicConfig(filename="./project.log",
        level=logging.INFO,
        format="\n[job.daemon][%(asctime)s][%(levelname)s]:\n%(message)s",
        datefmt='%y/%m/%d %H:%M:%S')

INFILEPATH = os.path.abspath("./infile")
OUTFILEPATH = os.path.abspath("./outfile")

class JobAtom():
    '''atom class of all jobs'''
    def __init__(self, pid, bundle):
        self.pid = pid
        if type(bundle.execute) is str:
            self.execute = bundle.execute
        elif type(bundle.execute) is list:
            self.execute = " ".join(bundle.execute)
        else:
            print "Jobs.execute should be a list or str!"

        self.is_cluster = bundle.is_cluster
        self.auto_run = bundle.auto_run
        self.keep_cpu_busy = bundle.keep_cpu_busy
        self.name = bundle.name
        self.input_str = bundle.to_string(pid)
        return

    def get_job_name(self):
        '''get the name of JobAtom object'''
        return "Job({}).{}".format(self.name, self.pid)

def construct_job_queue(to_do):
    '''construct JobAtom queue from Job class '''
    logging.info("Constructing the job queue...")
    job_queue = []
    pid = 0
    if os.path.exists(INFILEPATH):
        filelist = [int(elem.split('_')[-1]) for elem in os.listdir(INFILEPATH)]
        filelist.sort()
        if len(filelist)!=0:
            pid = filelist[-1]
    #bundle is class job
    for bundle in [elem for elem in to_do if elem.keep_cpu_busy == False]:
    #running the jobs doesn't use much cpu first
        for _ in range(0, bundle.duplicate):
            pid += 1
            job_queue.append(JobAtom(pid, bundle))

    for bundle in [elem for elem in to_do if elem.keep_cpu_busy == True]:
    #running the jobs doesn't use much cpu first
        for _ in range(0, bundle.duplicate):
            pid += 1
            job_queue.append(JobAtom(pid, bundle))
    logging.info("Constructed the job queue!")
    return job_queue

def check_status():
    ''' check the status of submitted jobs,
    if the job is done, remove it from PROCLIST so new job could be submitted'''
    for elemp in PROCLIST:
        if elemp[0].poll() is not None:
            PROCLIST.remove(elemp)
            logging.info(elemp[1].get_job_name()+" is ended!")
            print elemp[1].get_job_name()+" is ended..."
    return

def submit_job(job_atom):
    '''submit a job to cluster or your own computer'''

    if(os.path.exists(INFILEPATH) is not True):
        os.system("mkdir "+INFILEPATH)
    if(os.path.exists(OUTFILEPATH) is not True):
        os.system("mkdir "+OUTFILEPATH)

    homedir = os.getcwd()
    jobname = homedir.split("/")[-1]+"."+job_atom.name

    infile = INFILEPATH+"/_in_{}_{}".format(job_atom.name, job_atom.pid)
    outfile = OUTFILEPATH+"/out_{}_{}.txt".format(
        job_atom.name, job_atom.pid)
    jobfile = os.path.abspath("./_job_{}_{}.sh".format(
        job_atom.name, job_atom.pid))
    #write input file into ./infile folder
    f_job = open(infile,"w")
    f_job.write(job_atom.input_str)
    f_job.close()
    f_allinput = open(os.path.abspath("./all_input.log"),"a")
    f_allinput.write("Job ID: {}, Job name: {}\n".format(
            job_atom.pid, job_atom.name))
    f_allinput.write(job_atom.input_str)
    f_allinput.close()
    if job_atom.is_cluster:
        fjob = open(jobfile,"w")
        fjob.write("#!/bin/sh\n"+"#PBS -N "+jobname+"\n")
        fjob.write("#PBS -o "+homedir+"/Output\n")
        fjob.write("#PBS -e "+homedir+"/Error\n")
        fjob.write("cd "+homedir+"\n")
        fjob.write("echo "+infile+" | "+job_atom.execute)
        fjob.close()
        if job_atom.auto_run:
            os.system("qsub "+jobfile)
            os.system("rm "+jobfile)
            logging.info(job_atom.get_job_name+" submitted!")
        else:
            print "You have to run "+job_atom.get_job_name()+" by yourself!"
    else:
        if job_atom.auto_run:
            shellstr = "echo "+infile+" | "+job_atom.execute+" >> "+outfile
            proc = subprocess.Popen(shellstr, shell=True)
            if job_atom.keep_cpu_busy:
                PROCLIST.append((proc, job_atom))

            logging.info(job_atom.get_job_name()+" is started...")
            logging.info("input:\n"+job_atom.input_str)
            logging.info("PID:{}\n".format(proc.pid))
            print job_atom.get_job_name()+" is started..."
        else:
            print "You have to run "+job_atom.get_job_name()+" by yourself!"
    return

if __name__ == "__main__":
    logging.info("Jobs manage daemon is started...")
    JOBQUEUE = construct_job_queue(inlist.TO_DO)
    #print [e.keep_cpu_busy for e in JOBQUEUE]
    i = 0
    for ATOM in JOBQUEUE: 
        while ATOM.is_cluster is False and len(PROCLIST)>=inlist.CPU:
            check_status()
            time.sleep(inlist.SLEEP)
        submit_job(ATOM)

    check_status()
    while len(PROCLIST)!=0:
        time.sleep(inlist.SLEEP)
        check_status()

    logging.info("Jobs manage daemon is ended...")
