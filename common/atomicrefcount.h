// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/atomicrefcount.h
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

#if !defined(MPTOOLKIT_COMMON_ATOMICREFCOUNT_H)
#define MPTOOLKIT_COMMON_ATOMICREFCOUNT_H

#include <atomic>
#include "trace.h"

class AtomicRefCount
{
   public:
      AtomicRefCount();
      AtomicRefCount(int Initial);

      AtomicRefCount(AtomicRefCount&) = delete;
      AtomicRefCount& operator==(AtomicRefCount const&) = delete;

      // atomically increments the counter
      void operator++();

      // atomically decrements the counter.  Returns zero if the
      // new value of the counter is zero, otherwise returns a positive number.
      int operator--();

      // returns true if the counter is zero.
      bool is_zero() const;

      // returns true if the value of the counter is 1, false otherwise
      bool is_shared() const;

      int value() const { return Count.load(); }

   private:
      std::atomic<unsigned> Count;
};

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
   std::atomic_fetch_add_explicit(&Count, 1u, std::memory_order_relaxed);
}

inline
int AtomicRefCount::operator--()
{
   unsigned c = std::atomic_fetch_sub_explicit(&Count, 1u, std::memory_order_release);
   if (c == 1u)
   {
      std::atomic_thread_fence(std::memory_order_acquire);
      return 0;
   }
   return c-1;
}

inline
bool AtomicRefCount::is_zero() const
{
   return Count.load() == 0;
}

inline
bool AtomicRefCount::is_shared() const
{
   return Count.load() != 1;
}

class shared_counter
{
   public:
      shared_counter() { Count  = nullptr; }
      shared_counter(shared_counter const& Other) = default;
      shared_counter(shared_counter&& Other) : Count(Other.Count) { Other.Count = nullptr; }
      ~shared_counter() = default;
      shared_counter& operator=(shared_counter const& Other) = default;
      shared_counter& operator=(shared_counter&& Other)
      {
	 Count = Other.Count;
	 Other.Count = nullptr;
	 return *this;
      }

      void allocate(int InitialValue = 0)
      {
	 DEBUG_CHECK(!Count);
	 Count = new AtomicRefCount(InitialValue);
      }

      void deallocate()
      {
	 DEBUG_CHECK(Count);
	 delete Count;
	 Count = nullptr;
      }

      void operator++() const
      {
	 ++*Count;
      }

      int operator--() const
      {
	 return --*Count;
      }

      bool is_zero() const
      {
	 return Count->is_zero();
      }

      bool is_shared() const
      {
	 return Count->is_shared();
      }

   private:
      AtomicRefCount* Count;
};

#endif
