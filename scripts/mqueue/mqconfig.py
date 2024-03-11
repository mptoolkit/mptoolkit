#!/usr/bin/env python
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/mqueue/mqconfig.py
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

# local configuration for mq
# $Id$

import os

# ExecHost : the name of the frontend node that we ssh to run qsub
ExecHost='master'

# MQAddCommand : the name (including path, if necessary) of the mqadd command
MQAddCommand='mqadd'

# user name : obvious
UserName = os.getenv('USER')

# ScriptDir : this is where we store job scripts before they are executed
ScriptDir = os.path.expanduser('~/.mqueue/scripts')

# when we need to copy a job script, make a random name with this prefix
ScriptPrefix = 'mqscript'

# RunQueueFile : the file that stores the list of pending jobs
RunQueueFile = os.path.expanduser('~/.mqueue/runqueue')

# QueuedJobsFile : this stores the list of jobs that have been submitted to the queue,
# to work out what job script belongs to which job.  We could also use this to
# restart jobs that failed for some reason
QueuedJobsFile = os.path.expanduser('~/.mqueue/running')

# FailedJobsFile : if any jobs fail for some reason (ie. the qsub fails),
# the job is not kept in the queue but is instead added to this file.
# If the problem can be corrected then this file can be appended to the
# QueuedJobsFile
FailedJobsFile = os.path.expanduser('~/.mqueue/failed')

# QStatCommand : the full path to the qstat program
QStatCommand = '/opt/pbs/bin/qstat'

# QSubCommand : the full path we use to run qsub.  If we wanted to always
# use some extra options to qsub, we could add them here.
# This is only run from the frontend node, we don't need to ssh here.
QSubCommand = '/opt/pbs/bin/qsub'

# MaxQueueLength : maximum number of jobs we permit in the QUEUED state
MaxQueueLength = 1

# QueueDelay : after we successfully submit a job to the queue, wait this
# many seconds to see if it can run immediately so we can submit another
QueueDelay = 0.5

# ScriptStartupCommand : a string to add to the start of a job script.
# We use this to run the mq command on the frontend node, so that
# we can queue more jobs once one starts.
# We explicitly ignore the current job, meaning that, even if qstat
# sees it still as QUEUED, we don't count it and queue another.
ScriptStartupCommand = 'ssh -n '+ExecHost+' mq --ignore $PBS_JOBID'

# ScriptShutdownCommand : a string to add to the end of a job script.
# We use this to run the mq command on the frontend node, so that
# we can queue more jobs once one finishes.
ScriptShutdownCommand = 'ssh -n '+ExecHost+' mq --ignore $PBS_JOBID'
