// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/stackallocator.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
  stackallocator.h

  A fast stack-based allocator for tempoary memory allocations

  Created 2004-05-14 Ian McCulloch

  Enforces a strict LIFO ordering at runtime.

  TODO: This is not thread-safe.  Needs thread-specific data.
*/

#if !defined(STACKALLOCATOR_H_JHF473YR4783YHFP899PO)
#define STACKALLOCATOR_H_JHF473YR4783YHFP899PO

#include <cstddef>
#include "niftycounter.h"

namespace StackAlloc
{

// assume alignment of long is good enough for everything...
const std::size_t MinAlign = sizeof(long);

// free standing allocate/deallocate functions.
void* allocate(std::size_t size);
void deallocate(void* ptr, std::size_t size);

// helper to return a typed array **NOTE THIS DOES NO CONSTRUCTION OR DESTRUCTION**
template <typename T>
inline
T* allocate_type(std::size_t size)
{
   return static_cast<T*>(allocate(size*sizeof(T)));
}

template <typename T>
inline
void deallocate_type(T* ptr, std::size_t size)
{
   return deallocate(ptr, size*sizeof(T));
}

// internal use only
namespace detail
{

void StackAllocatorInit();
void StackAllocatorExit();

// declare a nifty counter to do the initialization and shutdown
namespace
{
  NiftyCounter::nifty_counter<StackAllocatorInit, StackAllocatorExit> PoolAllocCounter;
} // namespace

} // namespace detail

} // namespace StackAlloc

#endif
