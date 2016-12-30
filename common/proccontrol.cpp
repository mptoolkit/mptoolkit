// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/proccontrol.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
//
// The following symbol used to be needed is needed on Linux, so
// that sa_sigaction is defined according to POSIX.
//#ifdef _POSIX_C_SOURCE
//#undef _POSIX_C_SOURCE
//#endif
//#define _POSIX_C_SOURCE 199309L

#include "proccontrol.h"
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>

#include <cstdarg>
#include <iostream>
#include <iomanip>
#include "common/stringutil.h"

#if defined(CONFIG_PBS_WALLTIME)
#include "walltime.h"
#endif

#if defined(DEBUG_DGEMM)
#include "common/blas3f.h"
#include "blas3f.cpp"
#endif

namespace
{

// globals for async checkpointing

// these are set on a call to AsyncCheckpoint()
sig_atomic_t volatile AsyncCheckpointFlag = false;  // this may be set asynchronously in a signal handler
char AsyncCheckpointReason[1024] = "";
int AsyncCheckpointCode;
double AsyncCheckpointTime;

int TestAsyncCheckpointCalls = 0; // records the number of calls to TestAsyncCheckpoint()
double LastTestAsyncCheckpointTime = ProcControl::GetCPUTime();
double MaxTestAsyncCheckpointLag = 0; // maximum lag between calls to TestAsyncCheckpoint()
double AsyncCheckpointLagSumSquared = 0; // sum of the squared lag of async checkpoint.
                                         // The average lag is AsyncCheckpointLagSumSquared / total CPU time.

double MaxCPUTime = 0;
double MaxWallTime = 0;

double ProgStartWallTime = ProcControl::GetWallTime();

double PreviousCPUTime = 0;
double PreviousElapsedTime = 0;

   //bool BacktraceDone = false; // only do a single backtrace in a multithread program

#if defined(CONFIG_PBS_WALLTIME)
double ProgStartPBSTime = pbs_timeused();
double PBSWallTimeBound = 0;
double MaxPBSTime = 0;
#endif

// path name and executable name, for doing a backtrace
char* tracePath;
char* traceName;

// helper function to printf directly to a file descriptor.
// It is safe to call this from inside a signal handler.
// UPDATE 2005-09-02: it is actually still not safe to call this
// from a signal handler, as vsnprintf() is not async-signal safe.
// http://www.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_04.html#tag_02_04_03
// However in practice it seems much safer than other *printf
// functions, which very commonly call malloc() etc.
int fdprintf(int fd, char const* format, ...)
{
   std::va_list ap;
   char buf[1024];
   va_start(ap, format);
   int Size = vsnprintf(buf, 1024, format, ap);
   if (Size > 1024) Size = 1024;
   int n = write(fd, buf, Size);
   va_end(ap);
   return n;
}

// this is our signal handler for SIGTERM
extern "C"
void HandleTerm(int, siginfo_t* SigInfo, void* Context)
{
   ProcControl::AsyncCheckpointSignal(ProcControl::ReturnSIGTERM, "Caught SIGTERM");
}

// this is our signal handler for SIGINT
extern "C"
void HandleInt(int, siginfo_t* SigInfo, void* Context)
{

   ProcControl::AsyncCheckpointSignal(ProcControl::ReturnSIGINT, "Caught SIGINT");
}

// we can also catch SIGUSR1/SIGUSR2
extern "C"
void HandleUsr1(int, siginfo_t* SigInfo, void* Context)
{
   ProcControl::AsyncCheckpointSignal(ProcControl::ReturnSIGUSR1, "Caught SIGUSR1");
}

extern "C"
void HandleUsr2(int, siginfo_t* SigInfo, void* Context)
{
   ProcControl::AsyncCheckpointSignal(ProcControl::ReturnSIGUSR2, "Caught SIGUSR2");
}

char const* StrSigNum(int Sig)
{
   switch (Sig)
   {
      // These are POSIX1.1
      case SIGHUP    : return "SIGHUP";
      case SIGINT    : return "SIGINT";
      case SIGQUIT   : return "SIGQUIT";
      case SIGILL    : return "SIGILL";
      case SIGABRT   : return "SIGABRT";
      case SIGFPE    : return "SIGFPE";
      case SIGKILL   : return "SIGKILL";
      case SIGSEGV   : return "SIGSEGV";
      case SIGPIPE   : return "SIGPIPE";
      case SIGALRM   : return "SIGALRM";
      case SIGTERM   : return "SIGTERM";
      case SIGUSR1   : return "SIGUSR1";
      case SIGUSR2   : return "SIGUSR2";
      case SIGCHLD   : return "SIGCHLD";
      case SIGCONT   : return "SIGCONT";
      case SIGSTOP   : return "SIGSTOP";
      case SIGTSTP   : return "SIGTSTP";
      case SIGTTIN   : return "SIGTTIN";
      case SIGTTOU   : return "SIGTTOU";

      // These are SUSv2
      case SIGBUS    : return "SIGBUS";
         //      case SIGPOLL   : return "SIGPOLL";
      case SIGPROF   : return "SIGPROF";
      case SIGSYS    : return "SIGSYS";
      case SIGTRAP   : return "SIGTRAP";
      case SIGURG    : return "SIGURG";
      case SIGVTALRM : return "SIGVTARLM";
      case SIGXCPU   : return "SIGXCPU";
      case SIGXFSZ   : return "SIGXFSZ";

      // Other signals are arch-specific
#if 0
      case SIGWINCH  : return "SIGWINCH";
      case SIGINFO   : return "SIGINFO";
      case SIGEMT    : return "SIGEMT";
#endif

      default        : return "unknown signal";
   }
}

char const* StrSigCode(int Sig, int Code)
{
   switch (Code)
   {
      case SI_USER       : return "user (kill, sigsend or raise)";
      //      case SI_KERNEL     : return "the kernel";
      case SI_QUEUE      : return "sigqueue";
      case SI_TIMER      : return "timer expired";
      case SI_MESGQ      : return "mesq state changed";
      case SI_ASYNCIO    : return "Async IO completed";
      //      case SI_SIGIO      : return "queued SIGIO";
   }

   switch (Sig)
   {
      case SIGILL :
         switch (Code)
         {
            // SIGILL
            case ILL_ILLOPC    : return "illegal opcode";
            case ILL_ILLOPN    : return "illegal operand";
            case ILL_ILLADR    : return "illegal addressing mode";
            case ILL_ILLTRP    : return "illegal trap";
            case ILL_PRVOPC    : return "privileged opcode";
            case ILL_COPROC    : return "coprocessor error";
            case ILL_BADSTK    : return "internal stack error";
         }
         break;
      case SIGFPE :
         switch (Code)
         {
            // SIGFPE
            case FPE_INTDIV    : return "integer divide by zero";
            case FPE_INTOVF    : return "integer overflow";
            case FPE_FLTDIV    : return "floating point divide by zero";
            case FPE_FLTOVF    : return "floating point overflow";
            case FPE_FLTUND    : return "floating point underflow";
            case FPE_FLTRES    : return "floating point inexact result";
            case FPE_FLTINV    : return "floating point invalid operation";
            case FPE_FLTSUB    : return "subscript out of range";
         }
         break;
      case SIGSEGV :
         switch (Code)
         {
            // SIGSEGV
            case SEGV_MAPERR   : return "address not mapped to object";
            case SEGV_ACCERR   : return "invalid permissions for mapped object";
         }
         break;
      case SIGBUS :
         switch (Code)
         {
            // SIGBUS
            case BUS_ADRALN    : return "invalid address alignment";
            case BUS_ADRERR    : return "non-existant physical address";
            case BUS_OBJERR    : return "object specific hardware error";
         }
         break;

      case SIGTRAP :
         switch (Code)
         {
            // SIGTRAP
            case TRAP_BRKPT    : return "process breakpoint";
            case TRAP_TRACE    : return "process trace trap";
         }
         break;

      case SIGCHLD :
         switch (Code)
         {
            // SIGCHLD
            case CLD_EXITED    : return "child has exited";
            case CLD_KILLED    : return "child was killed";
            case CLD_DUMPED    : return "child terminated abnormally";
            case CLD_TRAPPED   : return "traced child has trapped";
            case CLD_STOPPED   : return "child has stopped";
            case CLD_CONTINUED : return "stopped child has continued";
         }
         break;


         /*
      case SIGPOLL :
         switch (Code)
         {
            // SIGPOLL
            case POLL_IN       : return "data input available";
            case POLL_OUT      : return "output buffers available";
            case POLL_MSG      : return "input message available";
            case POLL_ERR      : return "i/o error";
            case POLL_PRI      : return "high priority input available";
            case POLL_HUP      : return "device disconnected";
         }
         break;
         */
   }
   return "unknown code";
}

// this is a handler for SEGV, BUS, etc ec
extern "C"
void HandleSegv(int Sig, siginfo_t* SigInfo, void* Context)
{
   fdprintf(2, "Caught signal %d %s \nerrno=%d (%s)\nGenerated by: %s\n",
           Sig, StrSigNum(SigInfo->si_signo), int(SigInfo->si_errno), strerror(SigInfo->si_errno),
           StrSigCode(SigInfo->si_signo, SigInfo->si_code));
   switch (SigInfo->si_signo)
   {
      case SIGILL  :
      case SIGSEGV :
      case SIGBUS  :
      case SIGFPE  :
         fdprintf(2, "Offending address: %p\n", (void*) SigInfo->si_addr);
         break;
   }
   //#if defined(PROCCONTROL_BACKTRACE)
   //   if (!BacktraceDone)
   //   {
   //      BacktraceDone = true;
   ProcControl::GenerateBacktrace();
   //   }
   sleep(10);
   //#endif

#if defined(DEBUG_DGEMM)
   if (Sig == SIGFPE)
   {
      fdprintf(2, "SIGFLT raised, debug output of last call to DGEMM follows:\n");
      BLAS::DGS::DebugPrintDgemm();
   }
#endif

   abort();
   //   signal(Sig, SIG_DFL);
   //   raise(Sig);
}

} // namespace

namespace ProcControl
{

MessageLogger::Logger ProcControlLog("ProcControl", std::cerr, 10);

Checkpoint::~Checkpoint()
{
}

void AsyncCheckpoint(int Code, std::string const& Reason)
{
   if (AsyncCheckpointFlag)
   {
         notify_log(10, ProcControlLog) << "Ignoring multiple asynchronous checkpoint, code = "
                                        << Code << '\n';
      return;
   }
   AsyncCheckpointFlag = true;

   AsyncCheckpointTime = GetCPUTime();
   notify_log(10, ProcControlLog) << "Asynchronous checkpoint triggered, code = " << Code
                                  << ", reason = " << Reason << '\n';
   AsyncCheckpointCode = Code;
   strncpy(AsyncCheckpointReason, Reason.c_str(), 1023);
}

void AsyncCheckpointSignal(int Code, char const* Reason)
{
   if (AsyncCheckpointFlag) return;
   AsyncCheckpointFlag = true;
   AsyncCheckpointTime = GetCPUTime();
   AsyncCheckpointCode = Code;
   strncpy(AsyncCheckpointReason, Reason, 1023);
}

void TestAsyncCheckpoint()
{
   double Now = GetCPUTime();
   ++TestAsyncCheckpointCalls;
   double Lag = Now - LastTestAsyncCheckpointTime;
   LastTestAsyncCheckpointTime = Now;
   MaxTestAsyncCheckpointLag = std::max(Lag, MaxTestAsyncCheckpointLag);
   AsyncCheckpointLagSumSquared += Lag * Lag;

   if (MaxCPUTime > 0 && GetCPUTime() > MaxCPUTime)
   {
      notify_log(10, ProcControlLog) << "Checkpoint at CPU time limit exceeded by "
                                     << (GetCPUTime() - MaxCPUTime) << '\n';
      throw Checkpoint(ReturnCPUTimeExceeded, "CPU time limit exceeded.");
   }

   if (MaxWallTime > 0 && GetElapsedTime() > MaxWallTime)
   {
      notify_log(10, ProcControlLog) << "Checkpoint at walltime limit exceeded by "
                                     << (GetWallTime() - MaxWallTime) << '\n';
      throw Checkpoint(ReturnCPUTimeExceeded, "Walltime limit exceeded.");
   }

#if defined(CONFIG_PBS_WALLTIME)
   if (MaxPBSTime > 0 && GetPBSTimeRemaining() <= 0)
   {
      double Remain = GetPBSTimeRemaining();
      notify_log(10, ProcControlLog) << "Checkpoint at PBS walltime limit exceeded by "
                                     << (-Remain) << '\n';
      throw Checkpoint(ReturnCPUTimeExceeded, "PBS Walltime time limit exceeded.");
   }
#endif

   if (AsyncCheckpointFlag)
   {
      notify_log(10, ProcControlLog) << "At asynchronous checkpoint region, time lag is "
                                     << (GetCPUTime() - AsyncCheckpointTime) << '\n';
      throw Checkpoint(AsyncCheckpointCode, AsyncCheckpointReason);
   }
}

void ClearAsyncCheckpoint()
{
   AsyncCheckpointFlag = false;
}

double GetCPUTime()
{
   ::rusage usage;
   int Err = ::getrusage(RUSAGE_SELF, &usage);
   CHECK(Err == 0);
   return usage.ru_utime.tv_sec + usage.ru_stime.tv_sec
     + 1.0E-6 * (usage.ru_utime.tv_usec + usage.ru_stime.tv_usec);
}

double GetCumulativeCPUTime()
{
   return GetCPUTime() + PreviousCPUTime;
}

double GetCumulativeElapsedTime()
{
   return GetElapsedTime() + PreviousElapsedTime;
}

void SetCPUTimeLimit(double n)
{
   MaxCPUTime = n;
   notify_log(10, ProcControlLog) << "CPU time limit is "
                                  << n << (n == 0 ? " (infinite)" : " seconds") << '\n';
}

void SetWallTimeLimit(double n)
{
   MaxWallTime = n;
   notify_log(10, ProcControlLog) << "Walltime limit is "
                                  << n << (n == 0 ? " (infinite)" : " seconds") << '\n';
}

long GetResidentMemory()
{
   ::rusage usage;
   int Err = ::getrusage(RUSAGE_SELF, &usage);
   CHECK(Err == 0);
   //   return usage.ru_maxrss;
   return usage.ru_idrss;
}

double GetWallTime()
{
   struct ::timeval tp;
   ::gettimeofday(&tp, NULL);
    return tp.tv_sec + 1.0E-6 * tp.tv_usec;
}

double GetElapsedTime()
{
   return GetWallTime() - ProgStartWallTime;
}

#if defined(CONFIG_PBS_WALLTIME)
double GetPBSTime()
{
   return pbs_timeused();
}

double GetPBSElapsedTime()
{
   return pbs_timeused() - ProgStartPBSTime;
}

double GetPBSTimeRemaining()
{
   double Now = GetWallTime();
   if (Now < PBSWallTimeBound) return PBSWallTimeBound - Now;

   double Remain = MaxPBSTime - GetPBSElapsedTime();
   PBSWallTimeBound = Now + Remain;
   return Remain;
}

void SetPBSWallTimeLimit(double n)
{
   MaxPBSTime = n;
   notify_log(10, ProcControlLog) << "PBS Walltime limit is "
                                  << n << (n == 0 ? " (infinite)" : " seconds") << '\n';
}
#endif

void GenerateBacktrace()
{
   char pidstr[50];
   sprintf(pidstr, "%d", (int) getpid());
   char opt1[3] = "-c";
   char opt2[10] = "backtrace";
   char* args[] = {opt1, opt2, tracePath, pidstr, NULL};

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

int CreateUniqueFile(std::string& Name, bool Unlink)
{
   Name = Name + "XXXXXX";
   char Buf[300];
   strncpy(Buf, Name.c_str(), 300-1);
   int Res = mkstemp(Buf);
   Name = Buf;
   if (Unlink && Res != -1)
   {
      //unlink(Buf);  // it seems that unlink() doesn't work here, but remove() does?  why?
      remove(Buf); // I take it back.  None of these lines works as I expected.
   }
   return Res;
}

void Initialize(std::string const& ProgName,
                double CumulativeCPUTime_, double CumulativeWallTime_,
                bool CatchTERM, bool CatchSEGV)
{
   notify_log(20, ProcControlLog) << "Installing signal handlers...\n";
   tracePath = static_cast<char*>(malloc(ProgName.size()+1));
   strcpy(tracePath, ProgName.c_str());
   //   tracePath = strdup(ProgName.c_str());
   traceName = ::strrchr(tracePath, '/');
   if (traceName) ++traceName;
   else traceName = tracePath;

   PreviousCPUTime = CumulativeCPUTime_;
   PreviousElapsedTime = CumulativeWallTime_;

   struct sigaction Action;
   sigemptyset(&Action.sa_mask);
   Action.sa_flags = SA_SIGINFO | SA_RESTART;

   if (CatchTERM)
   {
      Action.sa_sigaction = HandleTerm;
      sigaction(SIGTERM, &Action, NULL);
      Action.sa_sigaction = HandleUsr1;
      sigaction(SIGUSR1, &Action, NULL);
      Action.sa_sigaction = HandleUsr2;
      sigaction(SIGUSR2, &Action, NULL);
      Action.sa_sigaction = HandleInt;
      Action.sa_flags = SA_SIGINFO | SA_RESTART | SA_RESETHAND;
             // RESETHAND: on a second interrupt, revert to default behaviour (terminate)
      sigaction(SIGINT, &Action, NULL);
      Action.sa_flags = SA_SIGINFO | SA_RESTART;
            // put flags back to default before the next call to sigaction
   }

   // setup the segv handler
   if (CatchSEGV)
   {
      Action.sa_sigaction = HandleSegv;
      int ret = sigaction(SIGSEGV, &Action, NULL);
      ret = sigaction(SIGBUS, &Action, NULL);
      ret = sigaction(SIGFPE, &Action, NULL);
      if (ret)
      {
         perror("ProcControl::Initialize()");
         exit(1);
      }
   }
   if (CatchTERM || CatchSEGV) notify_log(20, ProcControlLog) << "signal handlers installed.\n";
}

void Shutdown()
{
   notify_log(20, ProcControlLog) << "ProcControl is shutting down.\n";
   notify_log(10, ProcControlLog) << "CPU time used: " << GetCPUTime() << '\n';
   notify_log(10, ProcControlLog) << "Elapsed time: " << GetElapsedTime() << '\n';
#if defined(CONFIG_PBS_WALLTIME)
   notify_log(10, ProcControlLog) << "PBS walltime: " << GetPBSElapsedTime() << '\n';
#endif

   if (TestAsyncCheckpointCalls > 0)
   {
      double AverageLag = AsyncCheckpointLagSumSquared / (2.0 * LastTestAsyncCheckpointTime);
      MaxTestAsyncCheckpointLag = std::max(MaxTestAsyncCheckpointLag, AverageLag);
      notify_log(10, ProcControlLog) << "Asynchronous checkpoint time lag average/max is "
                                     << AverageLag << " / " << MaxTestAsyncCheckpointLag << '\n';
   }
}

} // namespace ProcControl
