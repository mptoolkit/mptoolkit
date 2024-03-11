#!/usr/bin/env python
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/spectral/basic_input_output.py
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

def test_configuration(confdict):
    print('*********** C O N F I G   F I L E   L I S T ******************')
    errorcode=0
    try:
        print('Lower Frequency boundary: '+confdict['wmin'])
    except:
        print('### ERROR: Statement wmin(float) is missing or mispelled or something like that!')
        errorcode=1
        
    try:
        print('Upper Frequency boundary: '+confdict['wmax'])
    except:
        print('Statement wmax(float) is missing or mispelled or something like that!')
        errorcode=1
        
    try:
        print('Initial Frequency stepwidth: '+confdict['deltaw'])
    except:
        print('### ERROR: Statement deltaw(float) is missing or mispelled or something like that!')
        errorcode=1
        
    try:
        print('Systemtype is set to: '+confdict['System'])
    except:
        print('### ERROR: Statement System(string)[Cluster,RZ] is missing or mispelled or something like that!')
        errorcode=1
        
    try:
        print('Operator whose spectral density is calculated: '+confdict['Operator'])
    except:
        print('### ERROR: Statement Operator(string) is missing or mispelled or something like that!')
        errorcode=1
        
    try:
        print('Highest spectral-density for cutoff at the edges: '+confdict['Scutoff'])
    except:
        print('### ERROR: Statement Scutoff(float) is missing or mispelled or something like that!')
        errorcode=1
        
    try:
        print('Number of maximal-jobs run parallel: '+confdict['MaxJobs'])
    except:
        print('### ERROR: Statement MaxJobs(integer) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('Lattice-File: '+confdict['lattice-file'])
    except:
        print('### ERROR: Statement lattice-file(string) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('Wavefunction-file: '+confdict['wavefunction-file'])
    except:
        print('### ERROR: Statement wavefunction-file(string) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('Path for results: '+confdict['resultpath'])
    except:
        print('### ERROR: Statement resultpath(string) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('DMRG-Configuration file: '+confdict['dmrg_config-file'])
    except:
        print('### ERROR: Statement dmrg_config-file(string) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('List of k-values to be calculated: '+confdict['k-list'])
    except:
        print('### ERROR: Statement k-list(list with separator ,) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('Broading of peaks is: '+confdict['broading'])
    except:
        print('### ERROR: Statement broading(float) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('The size of the unit-cell is: '+confdict['unit-size'])
    except:
        print('### ERROR: Statement unit-list(integer) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('The user name is: '+confdict['user'])
    except:
        print('### ERROR: Statement user(string) is missing or mispelled or something like that!')
        errorcode=1

    try:
        print('The convergence criteria is set to: '+confdict['conv_cutoff'])
    except:
        print('### ERROR: Statement conv_cutoff(float) is missing or mispelled or something like that!')
        errorcode=1

    print('--------------------------------------------------------------')
    return errorcode

def read_configuration(filepath):
    try:
        try:
            print('... read configuration file ...')
            file=open(filepath,'r')
            lines=file.readlines()
        except:
            print('Dumm gelaufen')
    finally:
        config_dict={}
        for element in lines:
            if(element.endswith('\n')):
                element=element[:-1]
            splitting=element.split('=')
            config_dict[splitting[0]]=splitting[1]
        if test_configuration(config_dict)==0:
            return config_dict
        else:
            return -1
