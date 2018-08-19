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
      AtomicRefCount() noexcept;
      AtomicRefCount(int Initial) noexcept;

      AtomicRefCount(AtomicRefCount&) = delete;
      AtomicRefCount& operator==(AtomicRefCount const&) = delete;

      // atomically increments the counter
      void operator++() noexcept;

      // atomically decrements the counter.  Returns zero if the
      // new value of the counter is zero, otherwise returns a positive number.
      int operator--() noexcept;

      // returns true if the counter is zero.
      bool is_zero() const noexcept;

      // returns true if the counter is exactly 1.  This operation is synchronized
      // with decrementing the counter.  Do we need to synchronize with incrementing
      // the counter too?  In that case, operator++ would need release semantics on the
      // counter.  I don't think so, because if the count is currently 1, and one
      // thread wants to increase the count while another wants to call is_unique() (in
      // preparation for modifying whatever the counter is protecting) then there must
      // be some external synchronization anyway.
      bool is_unique() const noexcept;

      // returns the value of the counter; useful for debug only as there
      // is no synchronization
      int value() const noexcept { return Count.load(); }

   private:
      std::atomic<unsigned> Count;
};

inline
AtomicRefCount::AtomicRefCount() noexcept : Count(0)
{
}

inline
AtomicRefCount::AtomicRefCount(int Initial) noexcept : Count(Initial)
{
}

inline
void AtomicRefCount::operator++() noexcept
{
   std::atomic_fetch_add_explicit(&Count, 1u, std::memory_order_relaxed);
}

inline
int AtomicRefCount::operator--() noexcept
{
   // decrement the count with release semantics: prior read or write operations
   // cannot be reordered to occur below the atomic fetch
   unsigned c = std::atomic_fetch_sub_explicit(&Count, 1u, std::memory_order_release);
   if (c == 1u)
   {
      // If the counter has hit zero, then we are going to call the destructor soon.
      // Don't reorder any operations that are called in the destructor above the
      // acquire barrier.  (note that the memory barrier also imposes an ordering that
      // no read operations can be shifted below it, but I don't think we use that constraint).
      std::atomic_thread_fence(std::memory_order_acquire);
      return 0;
   }
   return c-1;
}

inline
bool AtomicRefCount::is_zero() const noexcept
{
   return Count.load() == 0;
}

inline
bool AtomicRefCount::is_unique() const noexcept
{
   return Count.load(std::memory_order_acquire) == 1;
}

class shared_counter
{
   public:
      shared_counter() { Count  = nullptr; }
      shared_counter(shared_counter const& Other) = default;
      shared_counter(shared_counter&& Other) noexcept : Count(Other.Count) { Other.Count = nullptr; }
      ~shared_counter() = default;
      shared_counter& operator=(shared_counter const& Other) = default;
      shared_counter& operator=(shared_counter&& Other) noexcept
      {
	 Count = Other.Count;
	 Other.Count = nullptr;
	 return *this;
      }

      // equivalent to *this = shared_counter()
      void reset()
      {
         Count = nullptr;
      }

      void allocate(int InitialValue = 1)
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

      void operator++() const noexcept
      {
	 ++*Count;
      }

      int operator--() const noexcept
      {
	 return --*Count;
      }

      bool is_zero() const noexcept
      {
	 return Count->is_zero();
      }

      bool is_unique() const noexcept
      {
	 return Count->is_unique();
      }

   private:
      AtomicRefCount* Count;
};

#endif
