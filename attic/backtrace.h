// -*- C++ -*- $Id$

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
