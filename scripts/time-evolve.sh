#!/bin/bash
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# scripts/time-evolve.sh
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

#!/bin/sh

function die
{
   echo "fatal" $1
   exit 1
}

if [ $# -lt 3 ] ; then
   echo "usage: time-evolve <wavefunction> <time> <timestep> [other parameters ....]"
   exit 1
fi

wavefunc=$1
shift
timemax=$1
shift
timestep=$1
shift
time=0
cp $wavefunc $wavefunc.$time || die "wavefunction $wavefunc does not exist!"
if [ $timemax -lt 0 ] ; then
    while [ $time -gt $timemax ] ; do
            nexttime=`expr $time + $timestep`
            cp $wavefunc.$time $wavefunc.$nexttime
            mp-evolve-krylov -w $wavefunc.$nexttime -t $timestep $@ || die "problem running mp-evolve-krylov."
            time="$nexttime"
    done
else
    while [ $time -lt $timemax ] ; do
            nexttime=`expr $time + $timestep`
            cp $wavefunc.$time $wavefunc.$nexttime
            mp-evolve-krylov -w $wavefunc.$nexttime -t $timestep $@ || die "problem running mp-evolve-krylov."
            time="$nexttime"
    done
fi
