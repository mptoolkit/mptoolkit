// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/poolallocator.h
//
// Copyright (C) 2000-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  poolallocator.h

  Defines a template class for memory pool allocations.

  Created 2000-04-06 Ian McCulloch

  To use - pool_allocator<T>::allocate() returns a pointer to a memory block
  big enough for a T.  pool_allocator<T>::deallocate(void* ptr) frees the block.
  Class operator new/delete can just call allocate/deallocate directly.
  This is most easily done by IMPLEMENT_POOL_ALLOCATOR_NEW(class).

  The maximum size of an object that can be pool allocated is specified by
  MaxBlockAllocatorUnits.  This is in units of MinAlign, which is sizeof(long).
  The memory block granularity is DefaultAllocationUnits, again in units
  of MinAlign.

  This was renamed from blockallocator.h on 2000-06-07

  Two main modes of usage:
  template <class T> pool_allocator<T> is a class which contains
  two static functions to allocate() and deallocate() memory blocks
  suitable for storing an object of type T.
  The other mode requires the size to be a parameter of PoolAlloc::allocate() and PoolAlloc::deallocate().

  The maximum size object is MaxBlockAllocatorUnits * MinAlign bytes.

  Allocating zero-sized objects is legal; the actual allocation will be 1 unit.

  if NDEBUG is not defined, the allocation heap will be checked on program termination
  and leaked pointers will be reported.  If POOL_ALLOCATOR_VERBOSE is defined, the report contains
  more detailed information.
*/ 

#if !defined(MPTOOLKIT_COMMON_POOLALLOCATOR_H)
#define MPTOOLKIT_COMMON_POOLALLOCATOR_H

#include "niftycounter.h"
#include "trace.h"
#include <stddef.h>

#if defined(POOLALLOC_TRACE_DETAILED)
#define TRACE_POOLALLOC(Msg) TRACE(Msg)
#else
#define TRACE_POOLALLOC(Msg) DUMMY_TRACE(Msg)
#endif

namespace PoolAlloc
{

//typedef ct_assert<sizeof(void*) == sizeof(long)> PVoidMustHaveMaxAlignment;

const size_t MinAlign = sizeof(void*);

// maximum size of a fixed size allocator, in multiples of MinAlign
const size_t MaxBlockAllocatorUnits = 256;

// default allcation size in multiples of MinAlign.  This is expanded as necessary.
const int DefaultAllocationUnits = 10240;

// free standing allocate/deallocate functions.
void* allocate(size_t size);
void deallocate(void* ptr, size_t size);

// allocate/deallocate objects of type T
template <class T>
struct pool_allocator
{
   static T* allocate();
   static void deallocate(T* ptr);
};

#define IMPLEMENT_POOL_ALLOCATOR_NEW(Obj)					\
inline void* operator new(size_t size)						\
{										\
   DEBUG_CHECK(size == sizeof(Obj));						\
   void* Ptr = PoolAlloc::allocate(size);					\
   TRACE_POOLALLOC("PoolAlloc operator new")(typeid(Obj).name())(size)(Ptr);	\
   return Ptr;									\
}										\
										\
inline void operator delete(void* p)						\
{										\
   TRACE_POOLALLOC("PoolAlloc operator delete")(typeid(Obj).name())(p);		\
   PoolAlloc::deallocate(p, sizeof(Obj));					\
}

// For debugging, walks the heap and reports the count of allocated vs freed blocks
void walk_heap();

// internal use only
void PoolAllocatorInit();
void PoolAllocatorExit();

// declare a nifty counter to do the initialization and shutdown
namespace
{
  NiftyCounter::nifty_counter<PoolAllocatorInit, PoolAllocatorExit> PoolAllocCounter;
} // namespace

// inlines

template <class T>
inline
T* pool_allocator<T>::allocate()
{
   return static_cast<T>(PoolAlloc::allocate(sizeof(T)));
}

template <class T>
inline
void pool_allocator<T>::deallocate(T* ptr)
{
   PoolAlloc::deallocate(ptr, sizeof(T));
}

} // namespace PoolAlloc

#endif
