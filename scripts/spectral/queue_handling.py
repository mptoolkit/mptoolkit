#Author: Sebastian Smerat
import os,sys,time

def send_job_to_queue(configlist,k,w):
    os.system('qsub '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'_'+w+'.batch > ./qsub.log')
    file=open('./qsub.log','r')
    jobnumber=file.read(5)
    print("Job started with jobnumber: "+jobnumber)
    file.close()
    time.sleep(1)
    return jobnumber

def cluster_free(configlist):
    free_test=False
    os.system('qstat > qstat.log')
    file=open('./qstat.log','r')
    file.readline() # Skip the first
    file.readline() # two uninteresting lines of qstat
    count_queued=0
    for line in file:
        job=line.split(' ')
        for x in range(60):
            try:
                job.remove('')
            except:
                break
        if job[4]=='Q' and job[2]==configlist['user']:
            count_queued=count_queued+1
    if int(count_queued)>=int(configlist['MaxQueuedJobs']):
        print('Cluster is full, queued jobs: '+str(count_queued))
    return int(count_queued)<int(configlist['MaxQueuedJobs'])


def give_finished_jobs(running_jobs_list,configlist):
    os.system('qstat > qstat.log')
    file=open('./qstat.log','r')
    file.readline() # Skip the first
    file.readline() # two uninteresting lines of qstat
    finished_jobs_list=[]
    cluster_list=[]
    for line in file:
        job=line.split(' ')
        for x in range(60):
            try:
                job.remove('')
            except:
                break
        cluster_list.append(job[0][:5])
    for x in running_jobs_list:
        if x.split(',')[2] not in cluster_list:
            finished_jobs_list.append(x)
    print("Now finished: "+str(finished_jobs_list))
    return finished_jobs_list









