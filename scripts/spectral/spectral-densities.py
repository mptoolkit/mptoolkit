#!/usr/bin/env python
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/spectral/spectral-densities.py
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

#!/usr/bin/python

# Author: Sebastian Smerat

# Import the modules here

import sys,os,time
import basic_input_output as files
import initials as inits
import job_handling as jobs
import queue_handling as queue
import data_evaluation as eval
import list_editing as list


#Reload the modules

reload(files)
reload(inits)
reload(jobs)
reload(queue)
reload(list)
reload(eval)


# Get the name prefix and the config-file and set variables

try:
    prefix=inits.get_and_print_prefix(sys.argv[1])
except:
    print('No prefix given!')


try:
    configfile=sys.argv[2]
except:
    print('No configfile given!')

if inits.check_config_file(sys.argv[2])==-1:
    sys.exit()

overall_convergation=False
running_jobs_list=[]
finished_jobs_list=[]
run_counter=0

# Load the configuration file

configlist=files.read_configuration(configfile)
if configlist==-1:
    print('### ERROR during reading the configuration-file\n### Program will be stopped\n')
    sys.exit()


# Create the initial job grid

list_of_new_jobs = inits.create_init_joblist(configlist['k-list'],configlist['wmin'],configlist['wmax'],configlist['deltaw'])


# Check the lattice and wave file, create directories, the k_operator and the lanczos vectors

if inits.check_files(configlist)==-1:
    sys.exit()
inits.create_directory_structure(configlist)
inits.create_k_operator(configlist)
inits.create_lanczos_vectors(configlist)
wmax_new=inits.get_new_wmax(configlist)
configlist['wmax']=wmax_new
print(configlist['wmax'])
wmax_dict=inits.make_wmax_dict(configlist)
wmin_dict=inits.make_wmin_dict(configlist)
print(str(wmax_dict))
print(str(wmin_dict))

##################
# Here start the main loop which finalizes, when the variable overall_convergation is set to True
while overall_convergation==False:
##################
    
    run_counter=run_counter+1

    print('********  B e g i n   o f   t h e   l o o p  **************************')
    print('This is run number: '+str(run_counter))

# Start Jobs from list_of_new_jobs and update the running_jobs_list

    started_jobs_list=jobs.start_jobs(list_of_new_jobs,prefix,configlist)
    running_jobs_list.extend(started_jobs_list)
    print("Actually running jobs: "+str(running_jobs_list))


# Check which jobs have finished during the last loop and append those to finished_jobs_list and remove them from running_jobs_list

    new_finished_job_list=queue.give_finished_jobs(running_jobs_list,configlist)
    for x in new_finished_job_list:
        running_jobs_list.remove(x)


# Get the data from the newly finished jobs

    finished_jobs_list.extend(eval.get_spectral_density_list(prefix,configlist,new_finished_job_list))
    print('All finished jobs at this time: '+str(finished_jobs_list))


# Here it is checked, whether jobs have converged, or not. When not, new jobs will be created

    for w0 in finished_jobs_list:
        w0split=w0.split(',')
        print(w0+':')

# Does w0 lie at the lower boundary of the spectrum?
        
        if float(w0split[1])-float(0.0001) <= float(wmin_dict[w0split[0]]):   # NEVER COMPARE TWO FLOATS ON EQUALITY => CAUSES PROBLEMS NOTHING BUT PROBLEMS
            if float(w0split[2]) < float(configlist['Scutoff']):
                print('The lower boundary has converged!')
            else:
                wmin_new=float(wmin_dict[w0split[0]])-float(configlist['deltaw'])
                list_of_new_jobs.insert(0,w0split[0]+','+str(wmin_new))
                wmin_dict[w0split[0]]=str(wmin_new)
#
# Does w0 lie at the upper boundary of the spectrum?
#
        elif float(w0split[1])+float(0.0001) >= float(wmax_dict[w0split[0]]):
            print('wir sind in der Schleife')
            if float(w0split[2]) < float(configlist['Scutoff']):
                print('The upper boundary has converged!')
            else:
                print('wir sind in der nummer 2')
                wmax_new=float(wmax_dict[w0split[0]])+float(configlist['deltaw'])
                list_of_new_jobs.insert(0,w0split[0]+','+str(wmax_new))
                print('List of Jobs: '+ str(list_of_new_jobs))
                wmax_dict[w0split[0]]=str(wmax_new)
#
# If we are in the bulk of the spectrum, check for convergence, otherwise create new jobs
#
        else:
            wlower=list.get_next_smallest_w(finished_jobs_list,running_jobs_list,w0split[1],configlist,w0split[0])
            print('Next lower one: '+str(wlower))
            whigher=list.get_next_highest_w(finished_jobs_list,running_jobs_list,w0split[1],configlist,w0split[0])
            print('Next higher one: '+str(whigher))
            if wlower !=False and whigher != False:
                convergence=jobs.test_convergence(configlist,w0,wlower,whigher)
                if convergence==False:
                    new_jobs_for_conv=jobs.get_new_points(configlist,w0,wlower,whigher)
                    for x in new_jobs_for_conv:
                        if x not in list_of_new_jobs:
                            list_of_new_jobs.insert(0,x)
                            print(list_of_new_jobs)
#
#
# check for overall_convergation
 
    if list_of_new_jobs==[] and running_jobs_list==[]:
        overall_convergation=True
    
    eval.make_spectral_densities_files(prefix,configlist,finished_jobs_list)

    

##################
    print('------  E N D   O F   L O O P ----------------------------------------------')
    time.sleep(60)
# END OF LOOP
#################



