// -*- C++ -*- $Id$

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
