// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/backtrace.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  ShowBacktrace() function for generating a stack trace.

  Created 2004-07-27 Ian McCulloch

  This assumes that the 'backtrace' shell script is accessible
  from the current path.

  The backtrace requires the path to the executable.  On linux,
  this is simply the symlink /proc/self/exe.  On other operating
  systems, the executable path must be supplied.  Thus, the
  zero-argument version is only available on Linux.
*/

#if !defined(BACKTRACE_H_CSHUIREWHGIUHIU)

void ShowBacktrace(char const* ProgramName);

#if defined(__linux)
inline
void ShowBacktrace() { ShowBacktrace("/proc/self/exe"); }
#endif

#endif
