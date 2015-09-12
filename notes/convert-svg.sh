#!/bin/bash
for i in *.svg ; do 
   s=${i%.svg}.pdf
   if [ $s -ot $i ] ; then
      inkscape -f $i -D -A ${i%.svg}.pdf
   fi
done
