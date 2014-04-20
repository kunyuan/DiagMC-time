'''This is the input file of all jobs. 
   You have to add new job objects to TO_DO list
   if you want to run simulation.'''
import job_class as job
CPU = 4
SLEEP = 5    #check job status for every SLEEP seconds
TO_DO = []

#common dictionary for all jobs
com_dict={
    "Lx" :  4,
    "Ly" :  4,
    "Jcp" :  1.0,
    "Beta" :  0.9,
    "Order" :  4,
    }


# monte carlo job defintion
mc_dict={
    "__Execute" : "./gamma3.exe",
    "__Duplicate" : 3,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "IsLoad" : False,
    "Reweight" : [1,1,1,1],
    #"ReadFile" : "0.90_1_coll",
    "Sample" : 1000000,
    "Sweep" : 10,
    "Toss" : 1000,
    "Worm/Norm" : 0.5 
    }
mc_dict.update(com_dict)
TO_DO.append(job.JobMonteCarlo(mc_dict))

# self consist loop job definition
sc_dict={
    "__Execute" : ["python", "./run_loop.py"],
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : True, 
    "IsLoad" : True,
    "ReadFile" : "0.90_1_coll",
    }
sc_dict.update(com_dict)
TO_DO.append(job.JobConsistLoop(sc_dict))

# output loop job definition
ol_dict={
    "__Execute" : ["python", "./run_loop.py"],
    "__Duplicate" : 1,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "IsLoad" : True,
    "ReadFile" : "0.90_1_coll",
    }
ol_dict.update(com_dict)
TO_DO.append(job.JobOutputLoop(ol_dict))

# output numerical integration job definition
ni_dict={
    "__Execute" : ["./gamma3.exe"],
    "__Duplicate" : 0,
    "__IsCluster" : False,
    "__AutoRun" : True,
    "IsLoad" : False,
    }
ni_dict.update(com_dict)
TO_DO.append(job.JobIntegration(ni_dict))

if __name__ == "__main__":
    for e in TO_DO:
        print e.ToString(1)+"\n"

