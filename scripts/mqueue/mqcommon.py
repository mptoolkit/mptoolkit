#!/usr/bin/python3
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/mqueue/mqcommon.py
#
# Copyright (C) 2012 Ian McCulloch <ian@qusim.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Research publications making use of this software should include
# appropriate citations and acknowledgements as described in
# the file CITATIONS in the main source directory.
#----------------------------------------------------------------------------
# ENDHEADER

# $Id$
# common functions for the manual queuing system

from mqconfig import *
import sys,os,socket,fcntl,string,time,struct,re,getopt

# the allowed options for qsub
QsubOptions = "a:A:c:C:e:hj:k:l:m:M:N:o:p:q:r:S:u:v:VW:z"

def AddToRunQueue(argv):
    # read the options, so we can get 'top' if necessary
    try:
        optlist,args = getopt.getopt(argv, QsubOptions, 'top')
    except getopt.GetoptError:
        print "mqsub: error: error parsing options"
        usage()
        sys.exit(2)
    InsertAtTop = False
    optcopy = optlist[:]
    for m in optcopy:
        if m == ('--top',''):
            InsertAtTop = True
            optlist.remove(('--top',''))
            
    OptionsLine = ' '.join(map(lambda x: x[0]+' '+x[1], optlist))
    OptionsLine += ' ' + ' '.join(args)
    # open the run queue file
    if not os.path.exists(os.path.dirname(RunQueueFile)):
        os.path.makedirs(os.path.dirname(RunQueueFile))
    if InsertAtTop:
        RunQueue = open(RunQueueFile, "r+")
        fcntl.lockf(RunQueue, fcntl.LOCK_EX)
        FileLines = RunQueue.readlines()
        RunQueue.seek(0)
        RunQueue.truncate(0)
        RunQueue.write(OptionsLine+'\n')
        RunQueue.writelines(FileLines)
        RunQueue.close()
    else:
        RunQueue = open(RunQueueFile, "a")
        fcntl.lockf(RunQueue, fcntl.LOCK_EX)
        RunQueue.write(OptionsLine+'\n')
        RunQueue.close() 
    print 'added job to runqueue: '+OptionsLine
    return

def SubmitJobs(Count):
    RunQueue = open(RunQueueFile, "r+")
    fcntl.lockf(RunQueue, fcntl.LOCK_EX)
    
    rq = RunQueue.readlines()
    QueuedJobs = []
    NumQueued=0
    FailedJobs = []

    QueuedJobsF = open(QueuedJobsFile, "a")
    fcntl.lockf(QueuedJobsF, fcntl.LOCK_EX)

    Error=False
    while Count > 0 and len(rq) > 0 and not Error:
        #print "submitting: ",rq[0].strip()
        JobIdOutput = os.popen(QSubCommand+' '+rq[0], 'r')
        JobId = string.strip(JobIdOutput.readline())
        if JobId=='':
            print "SubmitJobs(): error running: "+QSubCommand+' '+rq[0]
            FailedJobs.append(rq[0])
            #Error=True
        else:
            print "Job is queued: "+JobId
            QueuedJobs.append((JobId,rq[0]))
            QueuedJobsF.write(JobId+' '+rq[0])
            Count=Count-1

        NumQueued=NumQueued+1
        rq.pop(0)

    RunQueue.seek(0)
    RunQueue.truncate()
    RunQueue.writelines(rq)
    RunQueue.close()
    QueuedJobsF.close()

    RunQueue = open(RunQueueFile, "r+")
    rq = RunQueue.readlines()
    RunQueue.close()

    if FailedJobs:
        FailedJobsF = open(FailedJobsFile, "a")
        FailedJobsF.writelines(FailedJobs)
        FailedJobsF.close()

    if Error:
        sys.exit(2)
    return NumQueued

def QueuedCount(IgnoreList = []):
    Output = os.popen(QStatCommand)
    QStatOut = Output.readlines();
    ToIgnore = '|'.join(IgnoreList)
    Count = 0
    for i in QStatOut:
        if not (len(ToIgnore) > 0 and re.match(ToIgnore, i)) and \
               re.match(r'\S+\s+\S+\s+'+UserName+r'\s+\S+\s+Q\s',i):
            Count = Count + 1
    return Count
    
def RunningCount(IgnoreList = []):
    Output = os.popen(QStatCommand)
    QStatOut = Output.readlines();
    ToIgnore = '|'.join(IgnoreList)
    Count = 0
    for i in QStatOut:
        if not (len(ToIgnore) > 0 and re.match(ToIgnore, i)) and \
               re.match(r'\S+\s+\S+\s+'+UserName+r'\s+\S+\s+R\s',i):
            Count = Count + 1
    return Count

def TrySubmitJobs(IgnoreList = []):
    # I'm embarrased about this loop, too late at night to fix
    Suc=True
    while Suc:
        NumQueued = QueuedCount(IgnoreList)
        ToQueue = MaxQueueLength-NumQueued
        if ToQueue > 0:
            Suc=SubmitJobs(ToQueue) > 0
            #Suc=False
            if Suc:
                time.sleep(QueueDelay)
        else:
            Suc=False
        
