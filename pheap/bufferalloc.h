// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/bufferalloc.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  Buffer allocation for the persistent heap.

  Created 2002-07-14.

  Buffer management is centralized here to allow for easy changes to
  the buffer allocator, eg we might want to use mmap() or similar, in which
  case the pages need to be aligned appropriately.

  There are no public constructors for BufferAllocator, instead a mapping of
  PageSize -> BufferAllocator objects is maintained internally (such that there
  is at most one BufferAllocator object for each distinct page size), and
  returned by the named constructor BufferAllocator::GetAllocator().

  2004-02-24: added allocate_file() and read_file() functions for mmap-based
  file I/O.  Allocate_file() allocates a buffer as MAP_SHARED, thus any changes
  to the buffer are reflected in the file, once the buffer is deallocate()d.
  The file must be opened for reading and writing. If the mmap fails, allocate_file() will fail.
  Read_file() attempts to read a buffer from the specified file using mmap MAP_PRIVATE,
  and returns a read-only buffer.  Attempting to write to the buffer is likely to result
  in a SIGSEGV.  If the mmap() fails, an ordinary read() is attempted instead.
  The offset in the file must be a multiple of the OS pagesize.  Thus it is sufficient
  (but not necessary) that the offset be a multiple of the allocator page_size().
*/

#if ~defined(MPTOOLKIT_PHEAP_BUFFERALLOC_H)
#define MPTOOLKIT_PHEAP_BUFFERALLOC_H

#if defined(HAVE_CONFIG_H)
#include "config.h"    // for LARGEFILE configuration
#endif
#include <sys/types.h> // for off_t and size_t
#include <mutex>
#include <unordered_set>

// for large file support
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

class BufferAllocator
{
   public:
      size_t page_size() const { return PageSize; }

      unsigned char* allocate();
      void deallocate(unsigned char const* Buffer);

      unsigned char* allocate_file(int fd, off_t offset);
      unsigned char const* read_file(int fd, off_t offset);

      size_t get_page_size() const { return PageSize; }

      static BufferAllocator* GetAllocator(size_t PageSize);

      // For convenience of possible implementations,
      // the page size must be a multiple of PageGranularity,
      // currently set at 8KB (page size on Alpha architecture).
      static size_t const PageGranularity = 8192;

   private:
      BufferAllocator(size_t PageSize_);

      ~BufferAllocator(); // not implemented
      BufferAllocator(BufferAllocator const&); // not implemented
      BufferAllocator& operator=(BufferAllocator const);

      std::mutex PageBufferMutex;
      size_t PageSize;
      void* FirstFreePageBuffer;

      // So that we know which method to use when deallocating buffers, we keep
      // a hash of mmap'ed buffers.
      using MappedBuffersType = std::unordered_set<unsigned char const*>;
      MappedBuffersType MappedBuffers;
};

#endif
