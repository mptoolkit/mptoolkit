// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/linux_pthread_descr.h
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

#include </usr/include/pthread.h>

struct _pthread_descr_struct;
struct sigjmp_buf;

typedef _pthread_descr_struct* pthread_descr;


struct _pthread_descr_struct {
  pthread_descr p_nextlive, p_prevlive;
                                /* Double chaining of active threads */
  pthread_descr p_nextwaiting;  /* Next element in the queue holding the thr */
  pthread_descr p_nextlock;	/* can be on a queue and waiting on a lock */
  pthread_t p_tid;              /* Thread identifier */
  int p_pid;                    /* PID of Unix process */
  int p_priority;               /* Thread priority (== 0 if not realtime) */
  struct _pthread_fastlock * p_lock; /* Spinlock for synchronized accesses */
  int p_signal;                 /* last signal received */
  sigjmp_buf * p_signal_jmp;    /* where to siglongjmp on a signal or NULL */
  sigjmp_buf * p_cancel_jmp;    /* where to siglongjmp on a cancel or NULL */
  char p_terminated;            /* true if terminated e.g. by pthread_exit */
  char p_detached;              /* true if detached */
  char p_exited;                /* true if the assoc. process terminated */
  void * p_retval;              /* placeholder for return value */
  int p_retcode;                /* placeholder for return code */
  pthread_descr p_joining;      /* thread joining on that thread or NULL */
};

inline
std::ostream& operator<<(std::ostream& out, _pthread_descr_struct const& Desc)
{
  out << "p_nextlive         " << (void*) Desc.p_nextlive
      << "\np_prevlive         " << (void*) Desc.p_prevlive
      << "\np_nextwaiting      " << (void*) Desc.p_nextwaiting
      << "\np_nextlock         " << (void*) Desc.p_nextlock
      << "\np_tid              " << Desc.p_tid
      << "\np_pid              " << Desc.p_pid
      << "\np_priority         " << Desc.p_priority
      << "\np_signal           " << Desc.p_signal
      << "\np_signal_jmp       " << (void*) Desc.p_signal_jmp
      << "\np_cancel_jmp       " << (void*) Desc.p_cancel_jmp
      << "\np_terminated       " << (bool) Desc.p_terminated
      << "\np_detached         " << (bool) Desc.p_detached
      << "\np_exited           " << (bool) Desc.p_exited
      << "\np_retval           " << Desc.p_retval
      << "\np_retcode          " << Desc.p_retcode
      << "\np_joining          " << (void*) Desc.p_joining
      << "\n";
  return out;
}
