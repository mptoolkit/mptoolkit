// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pagebase.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  Generic interface for page files.

  Created 2004-02-13 Ian McCulloch

  An abstract interface for PageFile and FileSystem
*/

#if !defined(PAGEBASE_HJY4398Y98YF87Y9Y8WHYFGFAGIUDSG4)
#define PAGEBASE_HJY4398Y98YF87Y9Y8WHYFGFAGIUDSG4

#include "common/messagelogger.h"
#include "common/inttype.h"
#include "pstream/pstream.h"
#include "bufferalloc.h"
#include <stddef.h>

namespace pheap
{
extern MessageLogger::Logger PHeapLog;
} // namespace pheap

namespace PHeapFileSystem
{

using inttype::uint32;
using inttype::uint64;

class GenericPageFile
{
   public:
      GenericPageFile(size_t PageSize_, std::string const& Name_)
	: PageSize(PageSize_), Name(Name_) { Alloc = BufferAllocator::GetAllocator(PageSize_); }

      virtual ~GenericPageFile() = 0 {}

      // shuts down the page file.
      virtual void shutdown() = 0;

      // shuts down the page file.
      virtual void persistent_shutdown(uint64 UserData) = 0;

      // writes a buffer that is not currently on disk.  At this time,
      // a LocalPage is allocated to the PageInfo, and PageInfo->PF is set to this.
      // Precondition: PageInfo->PF == NULL.  Buffer was allocated via allocate_buffer().
      // Postcondition: The buffer is written to disk, and the buffer itself is deallocated,
      // ie, after this call, the buffer is owned by the PageFile implementation and must not
      // be further accessed or deallocated by the caller.
      // Returns the page number.
      virtual size_t write(unsigned char const* Buffer) = 0;

      // Reads a previously written page from disk.
      // Precondition: PageFile->PF == this, the page has previously been written
      // Postcondition: Returns a buffer containing the page.  The buffer should be deallocated with
      // deallocate_buffer().
      virtual unsigned char* read(size_t Page) = 0;
      
      // deallocates a previously written page.
      virtual void deallocate(size_t Page) = 0;

      BufferAllocator* get_allocator() const { return Alloc; }

      // helper functions for allocating/deallocating buffers.
      unsigned char* allocate_buffer() const { return get_allocator()->allocate(); }
      void deallocate_buffer(unsigned char const* Buffer) const { get_allocator()->deallocate(Buffer); }

      // returns the name of the page file
      std::string const& name() const { return Name; }

      // functions to get some statistics
      size_t get_page_size() const { return pageSize; }
      virtual size_t num_free_pages() const;       // the number of unused pages in the page file
      virtual size_t num_allocated_pages() const;  // the file size
      virtual int num_pages_written() const;       // number of pages written to disk, in the lifetime of this PageFile
      virtual int num_pages_read() const;          // number of pages read from disk, in the lifetime of this PageFile

      virtual void set_checkpoint_limit_kb(unsigned long Size) = 0;
      virtual unsigned long get_checkpoint_limit_kb() const = 0;

   private:
      PageFile(PageFile const&);            // not implemented
      PageFile& operator=(PageFile const&); // not implemented

      size_t PageSize;
      std::string Name;
      BufferAllocator* Alloc;
};

} // namespace PHeapFileSystem

#endif
