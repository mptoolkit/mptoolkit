#!/bin/bash
declare -i FailCount=0
declare -i Count=0
for i in $@; do
   echo running $i...
   eval $i || { echo "test $i failed with status $?" ; FailCount=$(($FailCount + 1)); }
   Count=$Count+1
done
if [[ $FailCount == 0 ]] ; then
   echo "All tests successful (total " $Count")"
elif [[ $FailCount == 1 ]] ; then
   echo "1 test failed!"
else
   echo "$FailCount tests failed!"
fi
