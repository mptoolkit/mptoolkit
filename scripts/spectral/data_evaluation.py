#!/usr/bin/env python
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/spectral/data_evaluation.py
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

import sys,os

def get_spectral_density_list(prefix,configlist,new_finished_job_list):
    finished_jobs=[]
    for set in new_finished_job_list:
        kw=set.split(',')
        k=kw[0]
        w=kw[1]
        file=open(str(configlist['resultpath'])+'/'+str(configlist['Operator'])+'_k_'+str(k)+'/'+str(prefix)+'_cv_'+str(k)+'_'+str(w)+'.sweep','r')
        lastline=''
        while 1:
            x=file.readline()
            if x=='':
                break
            else:
                lastline=x
        dos=lastline.split(' ')[8]
        file.close()
        finished_jobs.append(k+','+w+','+dos)
    return finished_jobs


def sort_kwS_list(element1,element2):
    element1_split=element1.split(',')
    element2_split=element2.split(',')

    if float(element1_split[1])<float(element2_split[1]):
        return -1
    if float(element1_split[1])==float(element2_split[1]):
        return 0
    if float(element1_split[1])>float(element2_split[1]):
        return 1


def make_spectral_densities_files(prefix,configlist,finished_jobs_list):
    print('********  E V A L U A T I O N  ********************************')
    
    k_list=configlist['k-list'].split(',')
    print('k-list: '+str(k_list))

    for k in k_list:
        kwS_list=[]
        sorted_list=[]
        for x in finished_jobs_list:
            x_split=x.split(',')
            if x_split[0]==k:
                kwS_list.append(x)

        kwS_list.sort(sort_kwS_list)
        
        file=open(configlist['resultpath']+'/'+prefix+'_'+configlist['Operator']+'_'+k+'.result','w')
        for x in kwS_list:
            x_split=x.split(',')
            file.write(x_split[0]+' '+x_split[1]+' '+x_split[2]+'\n')
        file.close()
        
    
    print('---------------------------------------------------------------')
