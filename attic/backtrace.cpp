// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/backtrace.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "backtrace.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

void ShowBacktrace(char const* Name)
{
   char pidstr[50];
   sprintf(pidstr, "%d", (int) getpid());
   char NameBraindeadCopy[256];
   strncpy(NameBraindeadCopy, Name, 255);
   char* const args[] = {"-c", "backtrace", NameBraindeadCopy, pidstr, NULL};

   // don't use system(), since that causes the child to ignore SIGINT.
   // instead we manually fork() and exec()
   pid_t pid = fork();
   if (pid == 0)
   {
      execv("/bin/sh", args);
   }
   else
   {
      waitpid(pid, NULL, 0);
   }
}
