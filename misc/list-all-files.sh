#!/bin/bash
# list all source files.  Use in conjunction with list-used-files
if [ $# -ne 1 ]; then
   echo "usage: list-files <source-path>"
   exit 1
fi

find $1 -name \*.h -o -name \*.cc -o -name \*.cpp | sort
