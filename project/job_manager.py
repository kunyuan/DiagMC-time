#!/usr/bin/python
import random
import os
import sys
import time
import subprocess
import logging

#IsCluster=True
IsCluster=False
TurnOnSelfConsist=False
cpu=4
#sourcedir="."
execute="./gamma3.exe"
loop_execute=['python','./run_self_consistent.py']
homedir=os.getcwd()
proclist=[]
loop_proc=[]
logging.basicConfig(filename=homedir+"/project.log",level=logging.INFO,format="\n[job.daemon][%(asctime)s][%(levelname)s]:\n%(message)s",datefmt='%y/%m/%d %H:%M:%S')
logging.info("Jobs manage daemon is started...")

def para_init():
    para=dict()
    para["IsForever"]=[]
    para["Sample"]=[]
    para["Sweep"]=[]
    para["Toss"]=[]
    para["IsLoad"]=[]
    para["Type"]=[]
    para["Lx"]=[]
    para["Ly"]=[]
    para["Jcp"]=[]
    para["Beta"]=[]
    para["Order"]=[]
    para["Reweight"]=[]
    para["Worm/Norm"]=[]
    para["ReadFile"]=[]
    para["Duplicate"]=[]
    para["IsCluster"]=[]
    return para

def check_status():
    time.sleep(10)
    flag = open("stop.inp", "r").readline().lstrip().rstrip()
    if(flag!='0'):
        logging.info("Jobs manage daemon is ended...")
        sys.exit()
    #os.system("clear")
    for elemp in proclist:
        if elemp[0].poll()!=None:
            proclist.remove(elemp)
            logging.info("Job "+str(elemp[1])+" is ended!")
            print "Job "+str(elemp[1])+" is ended..."
    return

def submit_jobs(para,execute,homedir):
    infilepath=homedir+"/infile"
    outfilepath=homedir+"/outfile"
    topstr="top -p"

    if(os.path.exists(infilepath)!=True):
        os.system("mkdir "+infilepath)
    if(os.path.exists(outfilepath)!=True):
        os.system("mkdir "+outfilepath)
    filelist=[int(elem.split('_')[-1]) for elem in os.listdir(infilepath)]
    filelist.sort()
    if(len(filelist)!=0):
        lastnum=filelist[-1]
    else:
        lastnum=0

    #if you want to iterate some parameter, add more "for" here
    for belem in para["Beta"]:
        for jelem in para["Jcp"]:

            if TurnOnSelfConsist:
                start=-1
            else:
                start=0

            for j in range(start,int(para["Duplicate"][0])):
                if j==-1: # when turn on self consisitent loop
                    pid=0
                else:
                    lastnum+=1
                    pid=lastnum  #the unique id and random number seed for job
                    pass
                
                infile="_in_"+str(pid)
                outfile="out_"+str(pid)+".txt"
                jobfile="_job_"+str(pid)+".sh"
                item=[]
                item.append(para["Lx"][0])
                item.append(para["Ly"][0])
                item.append(para["Toss"][0])
                item.append(para["Sample"][0])
                item.append(para["IsForever"][0])
                item.append(para["Sweep"][0])
                item.append(jelem)  #jcp
                item.append(belem)  #beta
                item.append(para["Order"][0])
                item.append(str(-int(random.random()*2**30)))   #Seed
                if j==-1:  # when turn on self consisitent loop
                    item.append(str(1))
                    item.append(str(1))
                else:
                    item.append(para["Type"][0])
                    item.append(para["IsLoad"][0])

                item.append(str(pid))
                item.append(para["ReadFile"][0])
                stri=" ".join(item)
                for eve in para["Reweight"]:
                    stri=stri+"\n"+eve
                stri=stri+"\n"+para["Worm/Norm"][0]+"\n\n"
                if IsCluster:
                    f=open(infilepath+"/"+infile,"w")
                    f.write(stri)
                    f.close()
                    f=open(homedir+"/all_input.log","a")
                    f.write("Job ID:"+str(pid))
                    f.write(stri)
                    f.close()
                    f=open(jobfile,"w")
                    f.write("#!/bin/sh\n"+"#PBS -N "+homedir.split("/")[-1]+"\n")
                    f.write("#PBS -o "+homedir+"/Output\n")
                    f.write("#PBS -e "+homedir+"/Error\n")
                    f.write("cd "+homedir+"\n")
                    if j>=0:
                        f.write("./"+execute+" < "+infilepath+"/"+infile)
                        f.close()
                        os.system("qsub "+jobfile)
                        os.system("rm "+jobfile)
                    else:
                        f.write(" ".join(loop_execute))
                        f.close()
                        os.system("qsub "+jobfile)
                        os.system("rm "+jobfile)
                    
                else: 
                    #print execute,infile
                    while True:
                        #print "proc:"+str(len(proclist))
                        if(len(proclist)>=cpu):
                            check_status()
                        else:
                            f=open(infilepath+"/"+infile,"w")
                            f.write(stri)
                            f.close()
                            f=open(homedir+"/all_input.log","a")
                            f.write(stri)
                            f.close()
                            f=open(infilepath+"/"+infile,"r")
                            g=open(outfilepath+"/"+outfile,"a")
                            if j>=0:
                                p=subprocess.Popen(execute,stdin=f,stdout=g)
                                proclist.append((p,pid))
                            else:
                                p=subprocess.Popen(loop_execute,stdin=f,stdout=g)
                                loop_proc.append(p)
                            f.close()
                            g.close()
                            topstr=topstr+str(p.pid)+","
                            f=open("./mytop.sh","w")
                            f.write(topstr[:-1])
                            f.close()                            
                            os.system("chmod +x ./mytop.sh")
                            if j==-1:
                                f=open("./kill_loop.sh","w")
                                f.write("kill -9 "+str(p.pid))
                                f.close()
                                os.system("chmod +x ./kill_loop.sh")
                            
                            logging.info("Job "+str(pid)+" is started...")
                            logging.info("input:\n"+stri)
                            print "Job "+str(pid)+" is started..."
                            break

    return

if len(sys.argv)>1:
    if sys.argv[1]=='l' or sys.argv[1]=='-l':
        logging.info("Self consistent loop is automatically started!")
        print "Self consistent loop is automatically started!"
        TurnOnSelfConsist=True
else:
    logging.info("Please take care of Self consistent loop by yourself!")
    print "Please take care of Self consistent loop by yourself!"
    TurnOnSelfConsist=False


inlist=open(homedir+"/inlist","r")
last=inlist.readline()[-2:]
if(last!="\r\n"):
    last="\n"

inlist.close()
inlist=open(homedir+"/inlist","r")
#parse inlist file
para=para_init()
for eachline in inlist:
    eachline=eachline.lstrip(' ').rstrip(last)
    
    if eachline=="":
        continue
    if eachline[0]=='%' and para["Sample"]!=0:
    #if you want to check the parameter in inlist, put the codes here
        flag=False
        for kelem in para.keys():
            if len(para[kelem])==0:
                logging.error(kelem+" is missing in inlist file!")
                flag=True
        if flag:
            sys.exit()
        if len(para['Reweight']) != int(para['Order'][0]):
            logging.error("Reweight doesn't match Order!")
            sys.exit()
        low=para['IsCluster'][0].lower()
        if low=="false" or low==".false.":
            IsCluster=False
        elif low=="true" or low==".true.":
            IsCluster=True
        else:
            logging.error("I don't understand '"+eachline+"'!")
            sys.exit()

        submit_jobs(para,execute,homedir)
        para=para_init()
        continue

    if eachline[0]!='#':
        try:
            key=eachline.split(":")[0].lstrip(' ').rstrip(' ')
            value=eachline.split(":")[1].lstrip(' ').rstrip(' ')
        except:
            logging.error("I don't understand '"+eachline+"'!")
            continue
        # parse list parameter
        if para.has_key(key):
            try:
                para[key]=[elem.lstrip(' ').rstrip(' ') for elem in value.split(',')]
            except ValueError:
                logging.error("Could not covert data "+key+":"+value)
                sys.exit()
            except:
                logging.error("Error happens!")
        else:
            logging.error("What is "+key+"?")
            sys.exit()

while True:
    if len(proclist)!=0:
        check_status()
    else:
        #print TurnOnSelfConsist,loop_proc
        if TurnOnSelfConsist and len(loop_proc)!=0:
            loop_proc[0].kill()
            logging.info("Loop daemon is killed.")
            print "Loop daemon is killed."
        break

logging.info("Jobs manage daemon is ended...")
