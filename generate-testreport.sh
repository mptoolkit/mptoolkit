#!/bin/bash
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# generate-testreport.sh
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
