// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/atomic_pthread.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
/*
  pthreads version of atomic.h

  Created 2002-08-12 Ian McCulloch
*/

#include "common/mutex.h"

class AtomicRefCount
{
   public:
      AtomicRefCount();
      AtomicRefCount(int Initial);

      // atomically increments the counter
      void operator++();

      // atomically decrements the counter.  Returns zero if the
      // new value of the counter is zero, otherwise returns a positive number.
      int operator--();

      // returns true if the counter is zero.
      bool is_zero() const;

      int value() const { pthread::mutex::sentry L(Lock); return Count; }

   private:
      AtomicRefCount(AtomicRefCount const&); // not implemented
      AtomicRefCount& operator=(AtomicRefCount const&); // not implemented

      int Count;
      mutable pthread::mutex Lock;
};

// it is assumed that in the pthreads case, whatever does the locking
// will provide memory barriers for us.
inline void memory_barrier()
{
}

inline void write_memory_barrier()
{
}

inline void read_memory_barrier()
{
}

inline
AtomicRefCount::AtomicRefCount() : Count(0)
{
}

inline
AtomicRefCount::AtomicRefCount(int Initial) : Count(Initial)
{
}

inline
void AtomicRefCount::operator++()
{
   pthread::mutex::sentry MyLock(Lock);
   ++Count;
}

inline
int AtomicRefCount::operator--()
{
   pthread::mutex::sentry MyLock(Lock);
   --Count;
   return Count;
}

inline
bool AtomicRefCount::is_zero() const
{
   pthread::mutex::sentry MyLock(Lock);
   return Count == 0;
}
