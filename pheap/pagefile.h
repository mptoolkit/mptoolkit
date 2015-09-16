// -*- C++ -*- $Id$

/*
  pagefile.h

  Raw page file interface

  Created 2002-07-12, but based on pheapfsv3

  A PageFile represents a sequence of pages on disk.
  Pages are fixed-size buffers that must be allocated with allocate_buffer(),
  and can then be written to disk or read back.

  The first page is reserved for the implementation, and is used to store metadata associated
  with the PageFile itself.

*/

#if !defined(PAGEFILE_H_FUH3897HY789RFH89HWE893487Y3489)
#define PAGEFILE_H_FUH3897HY789RFH89HWE893487Y3489

#if defined(HAVE_CONFIG_H)
#include "config.h" // for LARGEFILE configuration
#endif
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

class PageFileImpl;

class PageFile
{
   public:
      PageFile();

      // pseudo-constructor, creates a new PageFile, overwrites the old filename if it exists.
      void create(size_t PageSize, std::string const& FileName, bool Unlink = false);

      // pseudo-constructor, reopens an existing page file, returns the saved 'UserData' 
      // parameter from persistent_shutdown(int UserData).
      uint64 open(std::string const& FileName, bool ReadOnly = false);

      // shuts down the page file.
      void shutdown(bool Remove = false);

      // shuts down the page file.
      void persistent_shutdown(uint64 UserData);

      std::string const& name() const;

      // returns the version number of the page file.
      int version() const;

      // writes a buffer that is not currently on disk.  At this time,
      // a LocalPage is allocated to the PageInfo, and PageInfo->PF is set to this.
      // Precondition: PageInfo->PF == NULL.
      // Postcondition: The buffer is written to disk, and the buffer itself is deallocated,
      // ie, after this call, the buffer is owned by the PageFile implementation.
      size_t write(unsigned char const* Buffer);

      // Reads a previously written page from disk.
      // Precondition: PageFile->PF == this, the page has previously been written
      // Postcondition: Returns a buffer allocated with get_allocator().
      unsigned char const* read(size_t Page);
      
      // deallocates a previously written page.
      void deallocate(size_t Page);

      // If Page is in the last this->num_free_pages() pages on disk, then
      // physically move it to earlier in the file, thus shrinking the maximum length
      // of the page file.  Returns the new location (which may be the same as the
      // old location, if no move occurred).
      size_t try_defragment(size_t Page);

      BufferAllocator* get_allocator() const;

      // helper functions for allocating/deallocating buffers.
      unsigned char* allocate_buffer() const { return get_allocator()->allocate(); }
      void deallocate_buffer(unsigned char const* Buffer) const { get_allocator()->deallocate(Buffer); }

      // functions to get some statistics
      size_t get_page_size() const;
      size_t num_free_pages() const;       // the number of unused pages in the page file
      size_t num_allocated_pages() const;  // the file size
      int num_pages_written() const;       // number of pages written to disk, in the lifetime of this PageFile
      int num_pages_read() const;          // number of pages read from disk, in the lifetime of this PageFile

      void set_checkpoint_limit_kb(unsigned long Size);
      unsigned long get_checkpoint_limit_kb() const;

      void Debug();

   private:
      PageFile(PageFile const&);            // not implemented
      PageFile& operator=(PageFile const&); // not implemented

      PageFileImpl* Impl;
};

} // namespace PHeapFileSystem

#endif
