// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pheap/bufferalloc.cpp
//
// Copyright (C) 2004-2021 Ian McCulloch <ian@qusim.net>
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

#include "bufferalloc.h"
#include "common/trace.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <memory>
#include <map>
#include <errno.h>
#include <set>

namespace
{
   std::map<size_t, BufferAllocator*> AllocatorList;
   std::mutex AllocatorListMutex;

} // namespace

BufferAllocator::BufferAllocator(size_t PageSize_)
   : PageSize(PageSize_), FirstFreePageBuffer(NULL)
{
}

std::set<void*> FreeBuffers;

unsigned char*
BufferAllocator::allocate()
{
   std::lock_guard<std::mutex> Lock(PageBufferMutex);
   if (FirstFreePageBuffer)
   {
      unsigned char* Ret = static_cast<unsigned char*>(FirstFreePageBuffer);
      FirstFreePageBuffer = *static_cast<void**>(FirstFreePageBuffer);
#if !defined(NDEBUG)
      // poison the buffer
      memset(Ret, 0xAC, PageSize);
#endif
      FreeBuffers.erase(Ret);
      return Ret;
   }

   unsigned char* Ret = static_cast<unsigned char*>(malloc(PageSize));
   if (!Ret) throw std::bad_alloc();
#if !defined(NDEBUG)
    // poison the buffer
    memset(Ret, 0xAC, PageSize);
#endif

   return Ret;
}

void BufferAllocator::deallocate(unsigned char const* Buf)
{
   std::lock_guard<std::mutex> Lock(PageBufferMutex);
   void* Buffer = static_cast<void*>(const_cast<unsigned char*>(Buf));

   CHECK(FreeBuffers.count(Buffer) == 0);

   // see if this buffer was allocated using mmap()
   if (MappedBuffers.count(Buf))
   {
      int Ret = munmap(Buffer, PageSize);  // this should never fail
      CHECK(Ret == 0);
      MappedBuffers.erase(Buf);
      return;
   }

   // otherwise it was allocated using malloc()

#if !defined(NDEBUG)
    // poison the buffer
    memset(Buffer, 0xCB, PageSize);
#endif
   *static_cast<void**>(Buffer) = FirstFreePageBuffer;
   FirstFreePageBuffer = Buffer;

   FreeBuffers.insert(Buffer);
}

unsigned char* BufferAllocator::allocate_file(int fd, off_t offset)
{
   std::lock_guard<std::mutex> Lock(PageBufferMutex);
   void* Buf = mmap(NULL, PageSize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, offset);
   if (Buf == MAP_FAILED)
   {
      PANIC("allocate_file(): mmap() failed: ")(strerror(errno));
   }
   unsigned char* CBuf = static_cast<unsigned char*>(Buf);
   MappedBuffers.insert(CBuf);
   return CBuf;
}

unsigned char const* BufferAllocator::read_file(int fd, off_t offset)
{
   void* Buf = mmap(NULL, PageSize, PROT_READ, MAP_PRIVATE, fd, offset);

   if (Buf != MAP_FAILED)
   {
      unsigned char* CBuf = static_cast<unsigned char*>(Buf);
      std::lock_guard<std::mutex> Lock(PageBufferMutex);
      MappedBuffers.insert(CBuf);
      return CBuf;
   }
   // else the map failed, fall back to read()

   Buf = this->allocate();
   ssize_t Ret = pread(fd, Buf, PageSize, offset);
   if (Ret == -1)
   {
      PANIC("read_file(): read() failed: ")(strerror(errno));
   }
   else if (Ret != ssize_t(PageSize))
   {
      PANIC("read_file(): read() returned less than a full page: ")(strerror(errno));
   }
   return static_cast<unsigned char*>(Buf);
}

BufferAllocator*
BufferAllocator::GetAllocator(size_t PageSize)
{
   std::lock_guard<std::mutex> Lock(AllocatorListMutex);
   // make sure that PageSize is a multiple of the PageGranularity
   PRECONDITION(PageSize % PageGranularity == 0);

   BufferAllocator*& Alloc = AllocatorList[PageSize];
   if (Alloc == NULL) Alloc = new BufferAllocator(PageSize);
   return Alloc;
}
