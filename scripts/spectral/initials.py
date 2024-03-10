#!/usr/bin/python3
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/spectral/initials.py
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

import os,errno

def check_config_file(filepath):
    check=0
    if os.access(filepath,os.F_OK)==False:
        print(filepath+": the config-file does not exist")
        check=-1

def create_init_joblist(k_bound_list,wmin,wmax,deltaw):
    print('... Create initial jobs list ...')
    w=float(wmin)
    stop=float(wmax)
    step=float(deltaw)
    k_split_list=(k_bound_list.split(','))
    joblist=[]
    for k in k_split_list:
        while w<stop:
            joblist.append(k+','+str(w))
            w=w+step
            if (w <= 0.0000001 and w > 0):
                w=0
            if (w >= -0.0000001 and w < 0):
                w=0
        w=float(wmin)
    return joblist

def get_new_wmax(configlist):
    w=float(configlist['wmin'])
    stop=float(configlist['wmax'])
    step=float(configlist['deltaw'])
    while w<stop:
        w=w+step
    return w-step

def get_and_print_prefix(Arg1):
    if Arg1=="":
        print("Something wrong with prefix. Start program like ./spectral-densities.py <<prefix>>")
        return -1
    else:
        print("\n*********** J O B   N A M E **********************************")
        print("Name of the Job is: "+Arg1)
        print("--------------------------------------------------------------")
        return Arg1
    
def create_k_operator(configlist):
    var='mp-make-k -l '+configlist['lattice-file']+' -i '+configlist['Operator']+' -u '+configlist['unit-size']
    print('************* Create the k-Operators: ************************')
    os.system(var)
    print('--------------------------------------------------------------')

def create_lanczos_vectors(configlist):
    k_split_list=(configlist['k-list'].split(','))
    print("**************************************************************")
    print("Creation of lanczos vectors ...")
    os.system('mp-attr '+configlist['wavefunction-file']+' Energy > ./gs_energy.dat')
    file=open('./gs_energy.dat','r')
    gs_energy=file.readline()
    file.close()
    for k in k_split_list:
        var='mp-apply '+configlist['lattice-file']+' '+configlist['Operator']+'_k\('+k+'\)'+' '+configlist['wavefunction-file']+' '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'.wave'
        var2='mp-attr '+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k+'.wave GroundstateEnergy='+gs_energy
        os.system(var)
        os.system(var2)
    print("... done")
    print("--------------------------------------------------------------")

def create_directory_structure(configlist):
    print('... creating directory structure ...')
    k_split_list=(configlist['k-list'].split(','))
    for k in k_split_list:
        if os.access(configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k,1)==False:
            print('create ...'+configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k)
            os.makedirs(configlist['resultpath']+'/'+configlist['Operator']+'_k_'+k)

def check_files(configlist):
    check=0
    if os.access(configlist['lattice-file'],os.F_OK)==False:
        print(configlist['lattice-file']+": the lattice-file does not exist")
        check=-1
    if os.access(configlist['wavefunction-file'],os.F_OK)==False:
        print(configlist['wavefunction-file']+": the wavefunction-file does not exist")
        check=-1
    return check


def make_wmax_dict(configlist):
    k_split_list=(configlist['k-list'].split(','))
    wmax_dict={}
    for k in k_split_list:
        wmax_dict[k]=configlist['wmax']
    return wmax_dict

def make_wmin_dict(configlist):
    k_split_list=(configlist['k-list'].split(','))
    wmin_dict={}
    for k in k_split_list:
        wmin_dict[k]=configlist['wmin']
    return wmin_dict
    
