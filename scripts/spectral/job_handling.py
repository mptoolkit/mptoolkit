#!/usr/bin/env python
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/spectral/job_handling.py
#
# Copyright (C) 2012 Ian McCulloch <ian@qusim.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Reseach publications making use of this software should include
# appropriate citations and acknowledgements as described in
# the file CITATIONS in the main source directory.
#----------------------------------------------------------------------------
# ENDHEADER

import os,math
import queue_handling as queue

reload(queue)

def write_Header(prefix,configlist,k,w):
    file=open(configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'_'+w+'.batch','w')
    file.write("#PBS -l nodes=1:ppn=1\n#PBS -N "+prefix+"_"+k+"\n#PBS -m be\n#PBS -l mem=3000mb")
    file.write("\n\n# This File calculates the groundstate of a klm-model\n\n")
    file.close()

def setup_correction_vector_job(prefix,configlist,k,w):
    write_Header(prefix,configlist,k,w)
    var1='mp-gmres-init -H '+configlist['lattice-file']+':H -w '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'.wave'+' -l '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'.wave'+' -c '+configlist['dmrg_config-file']+' -o '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'/'+prefix+'_cv_'+k+'_'+w+' -F '+ w +' -B '+ configlist['broading']
    var2='mp-gmres-resume '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'/'+prefix+'_cv_'+k+'_'+w
    file=open(configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'_'+w+'.batch','a')
    file.write(var1+'\n')
    file.write(var2)
    file.close()

def start_jobs(jobslist,prefix,configlist):
    count=0
    started_jobs_list=[]
    print('********* T R Y   T O   S T A R T   J O B S ...***************')
    while jobslist!=[] and queue.cluster_free(configlist)==True:
        kw=jobslist.pop(0)
        k=kw.split(',')[0]
        w=kw.split(',')[1]
        count=count+1
        setup_correction_vector_job(prefix,configlist,k,w)
        jobnumber = queue.send_job_to_queue(configlist,k,w)
        started_jobs_list.append(k+','+w+','+str(jobnumber))
        try:
            os.system('rm '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'_'+w+'.batch')
        except:
            print("An error occured when trying to delete a batch file")
    print('--------------------------------------------------------------')
    return started_jobs_list


def test_convergence(configlist,w0,wlower,whigher):
    w_w0=w0.split(',')[1]
    w_wlower=wlower.split(',')[1]
    w_whigher=whigher.split(',')[1]
    S_w0=w0.split(',')[2]
    S_wlower=wlower.split(',')[2]
    S_whigher=whigher.split(',')[2]
    
    area_low_high=0.5*(float(S_wlower)+float(S_whigher))*(float(w_whigher)-float(w_wlower))
    area_low_null=0.5*(float(S_wlower)+float(S_w0))*(float(w_w0)-float(w_wlower))
    area_null_high=0.5*(float(S_w0)+float(S_whigher))*(float(w_whigher)-float(w_w0))

    area=math.fabs(area_low_high-area_low_null-area_null_high)

    if area<float(configlist['conv_cutoff']):
        print("These three are converged")
        return True
    else:
        print("These three have not converged")
        return False


def get_new_points(configlist,w0,wlower,whigher):
    k=w0.split(',')[0]
    w_w0=w0.split(',')[1]
    w_wlower=wlower.split(',')[1]
    w_whigher=whigher.split(',')[1]

    w1_new = 0.5*(float(w_w0)+float(w_wlower))
    w2_new = 0.5*(float(w_w0)+float(w_whigher))
    print('Neue Jobs bei: '+k+','+str(w1_new),k+','+str(w2_new) )
    return [k+','+str(w1_new),k+','+str(w2_new)]
