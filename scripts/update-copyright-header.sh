#!/bin/bash

if [ $# -lt 1 ]; then
   echo "usage: update-copyright-header <mode> [files...]"
   echo "if no files are specified, then recursively search for all files with existing headers."
   echo "available modes:"
   echo "   --show   : show each file with the new header"
   echo "   --files  : don't update headers, just list the files to be examined and quit"
   echo "   --diff   : show the difference from the old file to the new file"
   echo "   --commit : rewrite the header for each file"
   exit
fi

mode="$1"
shift

top_header="// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/"

default_copyright="// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>"

main_header="//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version."

# write it this way to avoid this file matching
header_end="//"
header_end="$header_end ENDHEADER"

main_header="$main_header
$header_end"

files="$@"

if [ -z "$files" ]; then
   files=$(grep -rls "$header_end" . | grep -v 'update-copyright-header' | grep -v '~')
fi

if [ -z "$files" ]; then
   echo "no files."
   exit
fi

if [ x"$mode" == x"--files" ]; then
   echo "Files:"
   echo "$files"
   exit
fi

# for processing all C++ files
#files=$(grep -rls '// -\*- C++ -\*-' . | grep -v 'update-copyright-header' | grep -v '~')

for i in $files ; do

   filename=${i:2}
   echo "Processing $filename"

   copyright="$(grep '^// Copyright' $i)"

   if [ -z "$copyright" ]; then
      echo "Copyrights: default"
      copyright="$default_copyright"
   else
      echo "Copyrights:"
      echo "$copyright"
   fi

   text=$(
#      if [[ ( $i == *.h ) || ( $i == *.cc ) || ( $i == *.cpp ) ]]; then
#         echo "// -*- C++ -*-"
#      fi

      echo "$top_header"
      echo "//"
      echo "// $filename"
      echo "//"
      echo "$copyright"
      echo "$main_header"
      if grep -lq "$header_end" $i ; then
         sed '0,/^\/\/ ENDHEADER$/d' $i
      else
         grep -v '// -\*- C++ -\*-' $i
      fi
   )

   if [ x"$mode" == x"--commit" ] ; then
      echo "$text" > $i
   elif [ x"$mode" == x"--diff" ] ; then
      echo "$text" | diff $i -
   else
      echo "$text"
   fi

   echo

done
