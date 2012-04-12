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
