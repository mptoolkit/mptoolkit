#!/bin/sh
#
# Script to modify the output of makedepend, so that
# it strips the leading directory and duplicates the
# target as X.d X.o (the intention being that the
# dependency file depends on the same files that the .o file does)
# 
# The dependencies are written to standard output
#
makedepend -a -f- -Y $@ 2>/dev/null \
	| sed 's/^[^:]*\///;s/\(.*\)\.o:/\1.d \1.o:/'
