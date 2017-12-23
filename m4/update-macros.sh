#!/bin/bash

if [ ! -d "$HOME/git/autoconf-archive/m4" ]; then
   echo "autoconf archive not found!"
   exit 1
fi

for i in *.m4 ; do
   if [ -f "$HOME/git/autoconf-archive/m4/$i" ]; then
      echo -n "checking $i ..."
      if diff "$i" "$HOME/git/autoconf-archive/m4/" > /dev/null; then
         echo "$i is up to date."
      else
         echo "updating $i"
         cp "$HOME/git/autoconf-archive/m4/$i" ./
      fi
   else
      echo "$i is not in the autoconf archive."
   fi
done
