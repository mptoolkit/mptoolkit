#!/bin/bash
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# fix-timestamp.sh
#
# Copyright (C) 2015 Ian McCulloch <ian@qusim.net>
# Copyright (C) 2015 Seyed Saadatmand <s.saadatmand@uq.edu.au>
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
# fix-timestamp.sh: prevents useless rebuilds after a "git clone", since git doesn't track timestamps
# aclocal-generated aclocal.m4 depends on locally-installed
# '.m4' macro files, as well as on 'configure.ac'
touch aclocal.m4
touch stamp-h.in
sleep 1
# autoconf-generated configure depends on aclocal.m4 and on
# configure.ac
touch configure
# so does autoheader-generated config.h.in
touch config.h.in
# and all the automake-generated Makefile.in files
touch Makefile.in
