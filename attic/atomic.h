// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/atomic.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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
  Atomic.h

  Portable atomic operations.

  Created 2002-07-11 Ian McCulloch

*/

#if !defined(ATOMIC_H_FUHUI34789347Y58Y7Y7834YWY8943)
#define ATOMIC_H_FUHUI34789347Y58Y7Y7834YWY8943

#include "common/mutex.h"

/*

  A note on memory barriers.
  In general, a memory barrier is required when decrementing a reference count, irrespective of the
  value of the count.  Consider the following scenario:

  Thread 1                                            Thread 2

  Access shared object
  Decrement reference count, returns non-zero         Decrement reference count, returns zero
                                                      destroy the shared object

  Thread 1 requires a write memory barrier before decrementing the reference count, to make sure that
  when thread 2 detects that the reference count has gone to zero, it knows that thread 1
  has no pending access in the shared object (and can therefore destroy it safely).

  Thread 2 requires a read memory barrier after decrementing the count, to prevent any reads
  required by the object destructor from been prefetched before the count goes to zero.
  The barrier must be after the count is decremented, to prevent the following unlikely scenario:

  Thread 1                                            Thread 2
                                                      Barrier
                                                      (prefetch some data in the shared object)
  Access shared object
  (overwriting thread 2's prefetched data)
  Barrier
  Decrement reference count, returns non-zero
                                                      Decrement reference count, returns zero
                                                      destroy the shared object, utilizing the
                                                      prefetched data

  If access to the shared object is protected by a mutex, then the mutex provides the necessary memory
  barrier for thread 1.  Thread 2 can destroy the object without synchronizing with the mutex only
  with an explicit read memory barrier.  Note that the barrier is needed after the decrement only in the
  case that the reference count is zero.

  The current implementation does not assume that a mutex is used to access the shared object,
  and so unconditionally requires a memory barrier before decrementing the reference count.


  For copy-on-write semantics, the reference count behaviour is a bit different.  In particular, it
  can be assumed that when a thread is writing to the object, no other threads have access to it.
  But synchronization is required:

  Thread 1                                            Thread 2

  (assuming reference count == 1)
  Write to object
  Write barrier                                       Increment reference count
                                                      Read barrier

  that is, a write memory barrier is required before Thread 2 increments the reference count.
  Thread 2 needs a read memory barrier.  Thread 1's memory barrier implies the need for some
  external synchronization.  This is not provided by AtomicCOWCount.

  Thread 1                                            Thread 2

  reads counter
  assigns new counter
  starts copy                                         increments count
  subtracts original count



The user view of AtomicRefCount:

class AtomicRefCount
{
   public:
      AtomicRefCount();                   // initializes counter to zero
      AtomicRefCount(int Initial);

      // atomically increments the counter
      void operator++();

      // atomically decrements the counter.  Returns zero if the
      // new value of the counter is zero, otherwise returns a positive number.
      // If zero is returned, then this operation also acts as a memory barrier.
      int operator--();

      // returns true if the counter is zero.  Should this act as a memory barrier?
      bool is_zero() const;

   private:
      AtomicRefCount(AtomicRefCount const&); // not implemented
      AtomicRefCount& operator=(AtomicRefCount const&); // not implemented
};

*/

//
// template <typename T> atomic
//
// A wrapper for a generic atomic access.  This is designed around counters etc.
// The generic version uses mutex's.
//

template <typename T>
class atomic
{
   public:
      atomic() {}
      atomic(T const& t) : value(t) { }

      T read() const { pthread::mutex::sentry Lock(Mut); return value; }
      void operator=(T const& t) { pthread::mutex::sentry Lock(Mut); value = t; }

      T operator++() { pthread::mutex::sentry Lock(Mut); return ++value; }
      T operator++(int) { pthread::mutex::sentry Lock(Mut); return value++; }

      T operator--() { pthread::mutex::sentry Lock(Mut); return --value; }
      T operator--(int) { pthread::mutex::sentry Lock(Mut); return value--; }

      T operator+=(T const& t) { pthread::mutex::sentry Lock(Mut); return value += t; }
      T operator-=(T const& t) { pthread::mutex::sentry Lock(Mut); return value -= t; }

      T exchange(T const& t) { pthread::mutex::sentry Lock(Mut); T temp = value; value = t; return temp; }

   private:
      atomic(atomic const&);  // not implemented
      atomic& operator=(atomic const&); // not implemented

      mutable pthread::mutex Mut;
      T value;
};

#if defined(MULTITHREAD)

// multithread variations, there is specialist code for
// alpha, and a generic pthread version

#if defined(__osf__) && defined(__alpha__)
#include "atomic_alpha.h"
#else
#include "atomic_pthread.h"
#endif

#else

// single thread, the UP version of mutex.h should
// collapse the pthreads locking down to no-ops

#include "atomic_pthread.h"

#endif

#endif
