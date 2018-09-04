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

#include "cuda-setup.h"

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if !defined(HAVE_CUDA)
#error "CUDA is required to compile this program!"
#endif

#include "cuda.h"
#include "common/environment.h"
#include "gpu_buffer.h"

namespace cuda
{

// cuda-setup.h implementation

bool is_cuda_enabled()
{
   return true;
}

int num_cuda_devices()
{
   return cuda::num_devices();
}

std::string num_cuda_devices_str()
{
   return cuda::num_devices_str();
}

std::vector<std::string> get_cuda_device_names()
{
   int n = cuda::num_devices();
   std::vector<std::string> Result(n);
   for (int i = 0; i < n; ++i)
   {
      cuda::device_properties d = cuda::get_device_properties(i);
      Result[i] = d.name();
   }
   return Result;
}

int const dev = getenv_or_default("MP_CUDA_DEVICE", 0);
//int const num_dev = cuda::num_devices();

int mp_cuda_device()
{
   int d = getenv_or_default("MP_CUDA_DEVICE", 0);
   if (dev == -1 || dev >= cuda::num_devices())
      return -1;
   return d;
}

int setup_cuda()
{
   if (dev == -1 || dev >= cuda::num_devices())
      return -1;
   cuda::set_device(dev);
   return dev;
}

int setup_cuda_thread()
{
   if (dev == -1 || dev >= cuda::num_devices())
      return -1;
   cuda::set_device(dev);
   return dev;
}

// cuda.h implementation

namespace detail
{

// stream management

std::mutex StreamFreeListMutex;
std::list<cudaStream_t> StreamFreeList;

int NumAllocatedStreams = 0;

cudaStream_t AllocateStream()
{
#if defined(CUDA_SYNCHRONIZE)
   return 0;
#else
   std::lock_guard<std::mutex> lock(StreamFreeListMutex);
#if 1
   // we don't need to find a stream that is unused; they all will be by the
   // time we get here
   std::list<cudaStream_t>::iterator sIter = StreamFreeList.begin();
   if (StreamFreeList.empty())
#else
   // search for a stream that isn't in use, we search in reverse order
   // just to simplify the logic
   std::list<cudaStream_t>::iterator sIter = StreamFreeList.end();
   cudaError_t e = cudaErrorNotReady;
   while (e == cudaErrorNotReady && sIter != StreamFreeList.begin())
   {
      --sIter;
      cudaError_t e = cudaStreamQuery(*sIter);
   }
   if (e == cudaErrorNotReady) // did we fail to find an unused stream?
#endif
   {
      // create some more streams
      DEBUG_TRACE("construcing some more streams")(NumAllocatedStreams);
      for (int i = 0; i < 10; ++i)
      {
	 cudaStream_t s;
	 check_error(cudaStreamCreate(&s));
	 StreamFreeList.push_back(s);
	 ++NumAllocatedStreams;
      }
      // return the last one
      sIter = StreamFreeList.end();
      --sIter;
   }
   cudaStream_t Result = *sIter;
   StreamFreeList.erase(sIter);
   return Result;
#endif
}

void FreeStream(cudaStream_t stream_)
{
   std::lock_guard<std::mutex> lock(StreamFreeListMutex);
   StreamFreeList.push_front(stream_);
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

// allocator for temporary GPU memory
blas::arena gpu_temp_arena(new BlockAllocator(DefaultBlockMultiple, false));

} // namespace detail

std::string
num_devices_str()
{
   int Result;
   cudaError_t e = cudaGetDeviceCount(&Result);
   if (e == cudaSuccess)
      return std::to_string(Result);
   if (e == cudaErrorNoDevice)
      return "0 (no devices)";
   if (e == cudaErrorInsufficientDriver)
      return "0 (insufficient driver)";
   return std::string("0 (") + cudaGetErrorString(e) + ")";
}

} // namespace cuda
