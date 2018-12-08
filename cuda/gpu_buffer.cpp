// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/gpu_buffer.cpp
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

#include "gpu_buffer.h"

namespace cuda
{

// container for storing the remains of gpu_buffer's that still have pending operations
// in the associated stream.  The purpose of this container is to keep track of the
// associated memory allocation so that we can free the memory asyncronously after the
// gpu_buffer has been destroyed.
using GpuBufferDeleteRecType = std::tuple<cuda::stream, blas::arena, void*, std::size_t>;
using GpuBufferDeleteListType = std::list<GpuBufferDeleteRecType>;

std::mutex GpuBufferPendingDeleteMutex;
GpuBufferDeleteListType GpuBufferPendingDelete;

void TryFlushGpuBuffers()
{
   std::lock_guard<std::mutex> lock(GpuBufferPendingDeleteMutex);
   GpuBufferDeleteListType::iterator I = GpuBufferPendingDelete.begin();
   int CountRunning = 0;
   while (I != GpuBufferPendingDelete.end())
   {
      if (std::get<0>(*I).is_running())
      {
	 ++I;
         ++CountRunning;
      }
      else
      {
	 std::get<1>(*I).free(std::get<2>(*I), std::get<3>(*I));
	 I = GpuBufferPendingDelete.erase(I);
      }
   }
   DEBUG_TRACE(CountRunning);
}

void SyncFlushGpuBuffers()
{
   std::lock_guard<std::mutex> lock(GpuBufferPendingDeleteMutex);
   GpuBufferDeleteListType::iterator I = GpuBufferPendingDelete.begin();
   while (I != GpuBufferPendingDelete.end())
   {
      std::get<0>(*I).synchronize();
      std::get<1>(*I).free(std::get<2>(*I), std::get<3>(*I));
      I = GpuBufferPendingDelete.erase(I);
   }
}

void AddToPendingDelete(cuda::stream& Stream, blas::arena& Arena, void* Ptr, std::size_t ByteSize)
{
   TRACE_CUDA("AddToPendingDelete")(Stream.raw_stream());
   {
      std::lock_guard<std::mutex> lock(GpuBufferPendingDeleteMutex);
      GpuBufferPendingDelete.emplace_back(std::move(Stream), std::move(Arena), Ptr, ByteSize);
   }
   SyncFlushGpuBuffers();
}

} // namespace cuda
