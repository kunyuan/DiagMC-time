#/usr/bin/python
import random
import os
import sys
import time
import subprocess

#IsCluster=True
IsCluster=False
cpu=4
#sourcedir="."
#execute="./test"
execute="./gamma3.exe"
homedir=os.getcwd()
proclist=[]

def para_init():
    para=dict()
    para["Block"]=[]
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
    return para

def check_status():
    time.sleep(10)
    for elemp in proclist:
        if elemp[0].poll()==0:
            print str(elemp[1])+" is done"
            proclist.remove(elemp)
            log=open("./logfile.txt","a")
            log.write("#"+str(elemp[1])+" job is ended at "+time.strftime("%Y-%m-%d %A %X %Z",time.localtime())+"\n")
            log.close()
    return

def submit_jobs(para,i,execute,homedir):
    infilepath=homedir+"/infile"
    outfilepath=homedir+"/outfile"
    if(os.path.exists(infilepath)!=True):
        os.system("mkdir "+infilepath)
    if(os.path.exists(outfilepath)!=True):
        os.system("mkdir "+outfilepath)

    #if you want to iterate some parameter, add more "for" here
    for belem in para["Beta"]:
        for jelem in para["Jcp"]:
            if IsCluster:
                infile="_in_selfconsist"
                item=[]
                item.append(para["Lx"][0])
                item.append(para["Ly"][0])
                item.append(para["Toss"][0])
                item.append(para["Sample"][0])
                item.append(para["Block"][0])
                item.append(para["Sweep"][0])
                item.append(jelem)  #jcp
                item.append(belem)  #beta
                item.append(para["Order"][0])
                item.append(str(-1))   #Seed
                item.append(str(2))
                item.append(str(1))
                item.append(para["ReadFile"][0])
                stri=" ".join(item)
                for eve in para["Reweight"]:
                    stri=stri+"\n"+eve
                stri=stri+"\n"+para["Worm/Norm"][0]+"\n\n"
                f=open(infilepath+"/"+infile,"w")
                f.write(stri)
                f.close()

            for j in range(0,int(para["Duplicate"][0])):
                i+=1
                pid=int(random.random()*2**30)  #the unique id and random number seed for job
                infile="_in_"+str(pid)
                outfile="out_"+str(pid)+".txt"
                jobfile="_job_"+str(pid)+".sh"
                item=[]
                item.append(para["Lx"][0])
                item.append(para["Ly"][0])
                item.append(para["Toss"][0])
                item.append(para["Sample"][0])
                item.append(para["Block"][0])
                item.append(para["Sweep"][0])
                item.append(jelem)  #jcp
                item.append(belem)  #beta
                item.append(para["Order"][0])
                item.append(str(-pid))   #Seed
                item.append(para["Type"][0])
                item.append(para["IsLoad"][0])
                item.append(para["ReadFile"][0])
                stri=" ".join(item)
                for eve in para["Reweight"]:
                    stri=stri+"\n"+eve
                stri=stri+"\n"+para["Worm/Norm"][0]+"\n\n"
                if IsCluster:
                    f=open(infilepath+"/"+infile,"w")
                    f.write(stri)
                    f.close()
                    f=open(homedir+"/all_input.txt","a")
                    f.write(stri)
                    f.close()
                    f=open(jobfile,"w")
                    f.write("#!/bin/sh\n"+"#PBS -N "+homedir.split("/")[-1]+"\n")
                    f.write("#PBS -o "+homedir+"/Output\n")
                    f.write("#PBS -e "+homedir+"/Error\n")
                    f.write("cd "+homedir+"\n")
                    f.write("./"+execute+" < "+infilepath+"/"+infile)
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
                            f=open(homedir+"/all_input.txt","a")
                            f.write(stri)
                            f.close()
                            f=open(infilepath+"/"+infile,"r")
                            g=open(outfilepath+"/"+outfile,"a")
                            p=subprocess.Popen(execute,stdin=f,stdout=g)
                            f.close()
                            g.close()
                            proclist.append((p,pid))
                            #print i,j
                            print str(pid)+" is started"
                            log=open("./logfile.txt","a")
                            log.write("#"+str(pid)+" job is started at "+time.strftime("%Y-%m-%d %A %X %Z",time.localtime())+"\n")
                            log.close()

                            break

    return i

#compile the source automatically

#filelist=os.listdir(sourcedir)
#sourcename=[elem for elem in filelist if elem[-3:]=="f90"]
#sourcename.sort()
#sourcename=sourcename[-1]
#compilerstr="ifort "+sourcedir+"/"+sourcename+" -O3 -o "+homedir+"/"+execute
#os.system(compilerstr)

inlist=open(homedir+"/inlist","r")
num=0
#parse inlist file
para=para_init()
for eachline in inlist:
    eachline=eachline.lstrip(' ').rstrip('\n')
    
    if eachline=="":
        continue
    if eachline[0]=='%' and para["Block"]!=0:
    #if you want to check the parameter in inlist, put the codes here
        flag=False
        for kelem in para.keys():
            if len(para[kelem])==0:
                print kelem+" is missing in inlist file!"
                flag=True
        if flag:
            sys.exit()
        if len(para['Reweight']) != int(para['Order'][0]):
            print "Reweight doesn't match Order!"
            sys.exit()
        num=submit_jobs(para,num,execute,homedir)
        para=para_init()
        continue

    if eachline[0]!='#':
        try:
            key=eachline.split(":")[0].lstrip(' ').rstrip(' ')
            value=eachline.split(":")[1].lstrip(' ').rstrip(' ')
        except:
            print "I don't understand '"+eachline+"'!"
            continue
        # parse list parameter
        if para.has_key(key):
            try:
                para[key]=[elem.lstrip(' ').rstrip(' ') for elem in value.split(',')]
            except ValueError:
                print "Could not covert data "+key+":"+value
                sys.exit()
            except:
                print "Error happens!"
        else:
            print "What is "+key+"?"
            sys.exit()

while True:
    if len(proclist)!=0:
        check_status()
    else:
        break

