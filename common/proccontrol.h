// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/proccontrol.h
//
// Copyright (C) 1999-2019 Ian McCulloch <ian@qusim.net>
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
  proccontrol.h

  traps the SIGTERM signal and uses it to set a global
  state flag, readable by the ShouldHalt() function.

  Created 1999-08-17 Ian McCulloch

  2001-12-04 : Modified to use sigaction() rather than signal(),
  added a handler for SIGSEGV, SIGBUS etc that writes a verbose
  error message before calling abort().
*/

#if !defined(MPTOOLKIT_COMMON_PROCCONTROL_H)
#define MPTOOLKIT_COMMON_PROCCONTROL_H

#include "messagelogger.h"
#include <string>

namespace ProcControl
{

// Log object used to output various bits of info.  Log stream defaults to
// std::cout.
//
// Levels used (default 10):
// 10     Startup/checkpoint/shutdown messages
//
extern MessageLogger::Logger ProcControlLog;

// Initializes the ProcControl.  This causes TestAsyncCheckpoint()
// to throw a Checkpoint exception if the current CPU time is > MaxCPUTime,
// or, if CONFIG_PBS_WALLTIME is defined and the current walltime is > maxPBSTime,
// or, if CatchTERM is set and either SIGINT, SIGTERM,SIGUSR1 or SIGUSR2 signal is caught.
// Passing a time of zero means infinite (no limit).  The SIGINT handler is one-shot;
// a second SIGINT signal reverts to the default behaviour and terminates the process.
// This also initializes the OpenMP state.
void Initialize(std::string const& ProgName,
                double CumulativeCPUTime_ = 0, double CumulativeWallTime_ = 0,
                bool CatchTERM = false, bool CatchSEGV = true);

// set the maximum CPU time, at which point an async checkpoint is signalled.
void SetCPUTimeLimit(double n);

// set the maximum elapsed time, at which point an async checkpoint is signalled.
void SetWallTimeLimit(double n);

// doesn't do much of importance; writes some statistics to ProcControlLog.
// should this be automated somehow?
void Shutdown();

// returns the cpu time of this process (system + user)
double GetCPUTime();

// returns the wall-time.  This is the standard UNIX time of the number
// of seconds since 1970-01-01.  Resolution is nominally 1 microsecond.
double GetWallTime();

// returns the elapsed time, in seconds, since the program started.
double GetElapsedTime();

// returns the cumulative CPU time, the cpu time of this process + previous time
double GetCumulativeCPUTime();

// returns the cumulative elapsed time, the elapsed time of this process + previous time
double GetCumulativeElapsedTime();

#if defined(CONFIG_PBS_WALLTIME)
// returns the absolute PBS walltime usage for the current job.
// (slow - requires connecting to PBS server)
double GetPBSTime();

// returns the PBS walltime since the program started.
// (slow - requires connecting to PBS server)
double GetPBSElapsedTime();

// returns the amount of remaining walltime, in seconds (from the bound supplied to Initialize())
// Note that only the sign is guaranteed accurate, the absolute value may be inaccurate.
// Only connects to the PBS server if the elapsed time is larger than the bound.
double GetPBSTimeRemaining();

// set the maximum PBS walltime, at which point an async checkpoint is signalled.
void SetPBSWallTimeLimit(double n);
#endif

// returns the amount of physical memory in use.  This is broken, doesn't work at all on linux.
long GetResidentMemory();

// forks, and calls the backtrace script to generate a stack trace
// at the current execution point.
void GenerateBacktrace();

// creates a unique temporary file using mkstemp.  On entry,
// Name is the initial component of the path name for the file.  Returns the file descriptor
// of the newly opened file, and sets Name' = actual name (= initial + 6 random characters).
// Note that the XXXXXX markers used by mkstemp(3) do not need to be added by the caller,
// but if only a directory is specified as the initial path name, the trailing slash is
// necessary!  If Unlink is true, the file is unlinked before returning.
int CreateUniqueFile(std::string& Name, bool Unlink = true);

// this class can be used to trigger a checkpoint, with throw Checkpoint(Code, Why)
class Checkpoint
{
   public:
      Checkpoint(int Code_, std::string const& Why_) : Code(Code_), Why(Why_) {}
      virtual ~Checkpoint();

      int ReturnCode() const { return Code; }
      std::string const& Reason() const { return Why; }

   private:
      int Code;
      std::string Why;
};

//
// Return codes for standard checkpoint conditions.  These are integers
// so that they can be used as exit codes.
//

int const ReturnAbnormal = 1;                // signals an abnormal return, restart is not possible
int const ReturnSIGTERM = 2;                 // checkpoint in response to SIGTERM
int const ReturnCPUTimeExceeded = 3;         // checkpoint in response to CPU time exceeded
int const ReturnSynchronousCheckpoint = 4;   // generic synchronous checkpoint
int const ReturnDiskLimitExceeded = 5;       // checkpoint in response to disk limit exceeded
int const ReturnMemoryLimitExceeded = 6;     // checkpoint in response to memory limit exceeded
int const ReturnSIGUSR1 = 7;                 // checkpoint in response to SIGUSR1
int const ReturnSIGUSR2 = 8;                 // checkpoint in response to SIGUSR2
int const ReturnSIGINT = 9;                  // checkpoint in response to SIGINT

//
// AsyncCheckpoint
//
// Flags that a checkpoint should occur at the next call to TestAsyncCheckpoint()
//

void AsyncCheckpoint(int Code, std::string const& Reason);

// A version of AsyncCheckpoint that is safe to call from a signal handler.
// This suppresses some possible information messages that might otherwise be written to the log.
void AsyncCheckpointSignal(int Code, char const* Reason);

//
// TestAsyncCheckpoint
//
// Call this function at strategic points where it is safe to checkpoint.
// If AsyncCheckpoint() has previously been called, or some other condition it set
// (eg maximum CPU time is exceeded) then a Checkpoint exception is thrown.
//

void TestAsyncCheckpoint();

//
// ClearAsyncCheckpoint
//
// Typically, on an asynchronous checkpoint the program will terminate and never
// call TestAsyncCheckpoint() again.  However, it is possible to reset the state
// and process multiple asynchronous checkpoints.  Any asynchronous checkpoints
// that occurred between the last checkpoint and this call will be ignored.
//

void ClearAsyncCheckpoint();

} // namespace ProcControl

#endif
