// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cuda.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "cuda.h"

namespace cuda
{

namespace detail
{

// stream management

std::mutex StreamFreeListMutex;
std::list<cudaStream_t> StreamFreeList;

cudaStream_t AllocateStream()
{
   std::lock_guard<std::mutex> lock(StreamFreeListMutex);
   if (StreamFreeList.empty())
   {
      for (int i = 0; i < 10; ++i)
      {
	 cudaStream_t s;
	 check_error(cudaStreamCreate(&s));
	 StreamFreeList.push_back(s);
      }
   }
   cudaStream_t s = StreamFreeList.back();
   StreamFreeList.pop_back();
   return s;
}

void FreeStream(cudaStream_t stream_)
{
   std::lock_guard<std::mutex> lock(StreamFreeListMutex);
   StreamFreeList.push_back(stream_);
}

// event management

std::mutex EventFreeListMutex;
std::list<cudaEvent_t> EventFreeList;

cudaEvent_t AllocateEvent()
{
   std::lock_guard<std::mutex> lock(EventFreeListMutex);
   if (EventFreeList.empty())
   {
      for (int i = 0; i < 10; ++i)
      {
	 cudaEvent_t s;
	 check_error(cudaEventCreateWithFlags(&s, cudaEventDisableTiming | cudaEventBlockingSync));
	 EventFreeList.push_back(s);
      }
   }
   cudaEvent_t s = EventFreeList.back();
   EventFreeList.pop_back();
   return s;
}

void FreeEvent(cudaEvent_t event_)
{
   std::lock_guard<std::mutex> lock(EventFreeListMutex);
   EventFreeList.push_back(event_);
}

} // namespace detail

} // namespace cuda
