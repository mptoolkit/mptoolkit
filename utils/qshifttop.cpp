// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// utils/qshifttop.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include <iostream>
#include <iomanip>
#include <stack>
#include <vector>

extern "C"
{

#include <pbs_error.h>
#include <pbs_ifl.h>

} // extern "C"

using namespace std;

// Aborts on a PBS error
void Panic(int Conn)
{
  char* ErrorMsg;
  if ((ErrorMsg = pbs_geterrmsg(Conn)) != NULL)
    cerr << "qshifttop: " << ErrorMsg << '\n';
  else
    cerr << "qshifttop: Error " << pbs_errno << '\n';

  pbs_disconnect(Conn);
  exit(2);
}

int main(int argc, char** args)
{
   if (argc < 2)
   {
     cerr << "usage: qshifttop jobid...\n";
     return 1;
   }

   int Conn = pbs_connect(NULL);
   if (Conn < 0)
   {
     cerr << "qshifttop: could not connect to PBS server.\n";
     return 2;
   }

   vector<char*> Jobs(&args[1], &args[argc]);

   while (!Jobs.empty())
   {
      char* JobID = Jobs.back();
      Jobs.pop_back();

      int JobIDLen = strlen(JobID);

      // select all jobs in the 'Q' state with user of $USER
      if (!getenv("USER"))
      {
         cerr << "qshifttop: error: $USER is not set!\n";
         exit(1);
      }

      attropl SelectUser = {NULL, ATTR_A, "", getenv("USER"), EQ};
      attropl SelectQueued = {&SelectUser, ATTR_state, "", "Q", EQ};
      char** QueuedJobs = pbs_selectjob(Conn, &SelectQueued, NULL);
      if (pbs_errno) Panic(Conn);

      // Probe the list of queued jobs until we find JobID, and push them onto a stack
      stack<char*> JobStack;
      while (*QueuedJobs && strncmp(*QueuedJobs, JobID, JobIDLen) != 0)
      {
	 JobStack.push(*QueuedJobs++);
      }
      // *QueuedJobs now contains the full name of jobID, or NULL if it wasn't in the list

      if (*QueuedJobs == NULL)
      {
	 cerr << "qshifttop: Job " << JobID << " is not in QUEUED state.\n";
      }
      else
      {
	 // Walk the stack of queued jobs and bubble JobID (*QueuedJobs) to the top.
	 // This might fail at some point if one or more jobs have started running since
	 // we called pbs_selectjob().  This is no problem, so we don't print an error message if this fails.
	 while (!JobStack.empty())
	 {
	    pbs_orderjob(Conn, *QueuedJobs, JobStack.top(), NULL);
	    if (pbs_errno) break;
	    JobStack.pop();
	 }
      }
   }

   pbs_disconnect(Conn);
   return 0;
}
