#!/bin/sh
# fix-timestamp.sh: prevents useless rebuilds after a "git clone", since git doesn't track timestamps
# aclocal-generated aclocal.m4 depends on locally-installed
# '.m4' macro files, as well as on 'configure.ac'
touch aclocal.m4
sleep 1
# autoconf-generated configure depends on aclocal.m4 and on
# configure.ac
touch configure
# so does autoheader-generated config.h.in
touch config.h.in
# and all the automake-generated Makefile.in files
touch Makefile.in
