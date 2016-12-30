// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pagefile_st.cpp
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

/*
  pagefile_st.cpp

  Single-threaded implementation of the page file interface.
  Although this is single-threaded, it is thread safe and reenterable.  maybe.

  Created in antiquity, Ian McCulloch
*/

#include "pheap/pheaperror.h"
#include "common/proccontrol.h"
#include "pstream/pfilestream.h"
#include "common/inttype.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <algorithm>
#include <string.h>
#include <cstdio>     // for remove
#include <set>
#include <iterator>

/*
  the version 1 format of the per-file metadata is:

  Page 0 is reserved, and contains

  uint64 magic
  uint32 Metadata version number
  uint32 PageSize
  uint32 NumAllocatedPages
  uint64 UserData

  All in XDR format.

  After all of the pages (ie. at offset PageSize * NumallocatedPages), the list of
  unused pages is written.  The free list is defined in PageFileImpl as std::set<size_t>,
  but it is streamed as uint32.

  UserData is an unsigned 64 bit number that is supplied by the caller.  This number is
  returned by a subsequent call to open().  The UserData field is a hook for higher level
  libraries to store their own metadata.

  Version 2 moves the free list to the end of the metadata:

  uint32 FreeListSize
  vector<uint32> FreeList number of additional pages
  vector<uint32> FreeList first part

  Followed by, for each page on the FreeList additional pages, another
  vector<uint32> FreeList
*/

namespace PHeapFileSystem
{

using PStream::ipfilestream;
using PStream::opfilestream;

using inttype::uint32;
using inttype::uint64;

// a 64-bit 'random' number that identifies binary page files.
uint64 PageFileMagic = 0x1fe03dc25ba47986LL;
uint32 PageFileMetadataVersion = 2;
int PageFileMinHeaderSize = 8+4+4+4+8+4+4+4;

class PageFileImpl
{
   public:
      PageFileImpl();

      void create(size_t PageSize_, std::string const& FileName_, bool Unlink = false, bool AllowOverwrite = true);

      uint64 open(std::string const& FileName_, bool ReadOnly);

      void shutdown(bool Remove = false);

      void persistent_shutdown(uint64 UserData);

      size_t write(unsigned char const* Buffer);

      unsigned char const* read(size_t Page);

      void deallocate(size_t Page);

      size_t try_defragment(size_t Page);

      BufferAllocator* get_allocator() const { return Alloc; }

      std::string const& name() const { return FileName; }

      int version() const { return MetaVersion; }

      size_t get_page_size() const { return PageSize; }
      size_t num_allocated_pages() const { return NumAllocatedPages; }
      size_t num_free_pages() const { return FreeList.size(); }
      int num_pages_written() const { return PagesWritten; }
      int num_pages_read() const { return PagesRead; }

      void set_checkpoint_limit_kb(unsigned long Size);
      unsigned long get_checkpoint_limit_kb() const;

      void Debug();

   private:
      size_t AllocatePage();
      bool IsOnFreeList(size_t Page);

      BufferAllocator* Alloc;

      int FD;  // the file descriptor
      std::string FileName;

      size_t PageSize;
      size_t NumAllocatedPages;   // number of pages in the file, whether they are currently used or not

      int MetaVersion;            // the version

      pthread::mutex FreeListMutex;
      std::set<size_t> FreeList;

      // Additional pages used to store the free list; we can deallocate these
      // when we close the file.
      std::list<inttype::uint32> FreeListAdditionalPages;

      bool ReadOnly;

      int PagesRead, PagesWritten;

      unsigned long PageCheckpointLimit;
};

PageFileImpl::PageFileImpl()
  : Alloc(NULL), FD(-1), PageSize(0), NumAllocatedPages(0), MetaVersion(PageFileMetadataVersion), ReadOnly(false),
    PagesRead(0), PagesWritten(0),
    PageCheckpointLimit(0)
{
}

void PageFileImpl::create(size_t PageSize_, std::string const& FileName_, bool Unlink, bool AllowOverwrite)
{
   CHECK(FD == -1);
   FileName = FileName_;
   PageSize = PageSize_;
   Alloc = BufferAllocator::GetAllocator(PageSize_);
   PagesRead = 0;
   PagesWritten = 0;
   NumAllocatedPages = 1;  // first page is reserved for implementation
   ReadOnly = false;
   FreeList.clear();
   MetaVersion = PageFileMetadataVersion;  // reset the version number to the default

   notify_log(20, pheap::PHeapLog) << "creating file " << FileName << '\n';
   int Flags = O_RDWR | O_CREAT | O_TRUNC;
   if (!AllowOverwrite)
      Flags |= O_EXCL;
   FD = ::open(FileName.c_str(), Flags, 0666);
   if (FD == -1)
   {
      throw pheap::PHeapCannotCreateFile(FileName, strerror(errno));
   }
   if (Unlink)
      ::unlink(FileName.c_str());
}

uint64 PageFileImpl::open(std::string const& FileName_, bool ReadOnly_)
{
   CHECK(FD == -1);
   FileName = FileName_;
   ReadOnly = ReadOnly_;
   PagesRead = 0;
   PagesWritten = 0;

   int OpenFlags = (ReadOnly ? O_RDONLY : O_RDWR);

   FD = ::open(FileName.c_str(), OpenFlags);
   if (FD == -1)
   {
      throw pheap::PHeapCannotOpenFile(FileName, strerror(errno));
   }

   // create a ipfilestream to read the metadata
   ipfilestream MetaIn(PStream::format::XDR);
   MetaIn.set_fd(FD);
   uint64 Magic = MetaIn.read<uint64>();

   if (Magic != PageFileMagic)
   {
      throw pheap::PHeapFileError(FileName_, "file format incorrect, cannot find magic number.");
   }

   uint32 Version = MetaIn.read<uint32>();
   notify_log(40, pheap::PHeapLog) << "PageFile " << FileName_ << " version number is " << Version << '\n';

   MetaVersion = Version;

   if (Version > 2)
   {
      throw pheap::PHeapFileError(FileName_, "file version mismatch, version is " +
                                  boost::lexical_cast<std::string>(Version) + " but expected version <= 2\n"
                                  "Probably this file was created with a newer version of the software.");
   }

   PageSize = MetaIn.read<uint32>();
   NumAllocatedPages = MetaIn.read<uint32>();
   uint64 UserData = MetaIn.read<uint64>();

   // Now we know the page size, get an allocator
   Alloc = BufferAllocator::GetAllocator(PageSize);

   if (Version == 1)
   {
      // Version 1 free list

      // This is a vector of uint32 stored at the end of the file.  We need
      // to convert from the default size_t
      off_t Offset = off_t(NumAllocatedPages) * off_t(PageSize);
      MetaIn.lseek(Offset, SEEK_SET);

      std::vector<inttype::uint32> TempFreeList;
      MetaIn >> TempFreeList;

      FreeList = std::set<size_t>(TempFreeList.begin(), TempFreeList.end());
   }
   else if (Version == 2)
   {
      //      TRACE("Found a version 2 file");
      // The version 2 free list stores it in allocated pages
      MetaIn >> FreeListAdditionalPages;
      std::vector<inttype::uint32> TempFreeList;
      MetaIn >> TempFreeList;
      FreeList = std::set<size_t>(TempFreeList.begin(), TempFreeList.end());

      for (std::list<inttype::uint32>::const_iterator I = FreeListAdditionalPages.begin();
           I != FreeListAdditionalPages.end(); ++I)
      {
         MetaIn.lseek(off_t(PageSize) * (*I), SEEK_SET);
         MetaIn >> TempFreeList;
         FreeList.insert(TempFreeList.begin(), TempFreeList.end());
      }
   }
   else
   {
      PANIC("Unsupported version number")(Version);
   }

   return UserData;
}

void PageFileImpl::shutdown(bool Remove)
{
   close(FD);
   FD = -1;
   if (Remove && !ReadOnly)
      std::remove(FileName.c_str());
}

// this is used by Debug() method
std::ostream& operator<<(std::ostream& out, std::set<std::size_t> const& l)
{
   std::copy(l.begin(), l.end(), std::ostream_iterator<std::size_t>(out, ", "));
   return out;
}

void PageFileImpl::persistent_shutdown(uint64 UserData)
{
#if 0
   // old version 1 format

   // write the free list at the end of the page file
   off_t Offset = off_t(NumAllocatedPages) * off_t(PageSize);
   lseek(FD, Offset, SEEK_SET);

   opfilestream MetaOut(PStream::format::XDR);
   MetaOut.set_fd(FD);

   // convert to uint32
   std::vector<inttype::uint32> TempFreeList(FreeList.begin(), FreeList.end());
   MetaOut << TempFreeList;

   // truncate the file at this point
   MetaOut.truncate();

   // The page size and number of allocated pages goes at the front
   MetaOut.lseek(0, SEEK_SET);
   MetaOut << PageFileMagic
           << PageFileMetadataVersion
           << uint32(PageSize)
           << uint32(NumAllocatedPages)
           << UserData;

   MetaOut.flush();
   MetaOut.close();
   FD = -1;


#else
   // Current version 2/3 format

   // We can free the old free list pages
   while (!FreeListAdditionalPages.empty())
   {
      this->deallocate(FreeListAdditionalPages.front());
      FreeListAdditionalPages.pop_front();
   }

   // Determine how many additional pages we need to store the free list
   uint32 FreeListBytes = FreeList.size() * 4;
   uint32 FirstPageBytes = PageSize - PageFileMinHeaderSize;
   while (FreeListBytes > FirstPageBytes)
   {
      // In this case we need to allocate additional pages
      FreeListAdditionalPages.push_back(this->AllocatePage());
      FreeListBytes = FreeList.size() * 4;
      FreeListBytes -= (PageSize - 4) * FreeListAdditionalPages.size();
      FirstPageBytes -= 4;   // for the 4 bytes to store the additional page number
   }

   opfilestream MetaOut(PStream::format::XDR);
   MetaOut.set_fd(FD);

   // truncate the file after the last allocated page
   off_t Offset = off_t(NumAllocatedPages) * off_t(PageSize);
   MetaOut.lseek(Offset, SEEK_SET);
   MetaOut.truncate();

   MetaOut.lseek(0, SEEK_SET);
   MetaOut << PageFileMagic
           << PageFileMetadataVersion
           << uint32(PageSize)
           << uint32(NumAllocatedPages)
           << UserData
           << FreeListAdditionalPages;

   std::vector<uint32> TempFreeList(FreeList.begin(), FreeList.end());
   uint32 Index = FirstPageBytes / 4;
   if (Index > TempFreeList.size())
      Index = TempFreeList.size();
   std::vector<uint32> FreeListSection(TempFreeList.begin(), TempFreeList.begin()+Index);
   MetaOut << FreeListSection;
   for (std::list<inttype::uint32>::const_iterator I = FreeListAdditionalPages.begin();
        I != FreeListAdditionalPages.end(); ++I)
   {
      MetaOut.lseek(off_t(PageSize) * (*I), SEEK_SET);
      uint32 NextIndex = Index + PageSize/4 - 1;
      if (NextIndex > TempFreeList.size())
         NextIndex = TempFreeList.size();
      FreeListSection = std::vector<uint32>(TempFreeList.begin()+Index,
                                            TempFreeList.begin()+NextIndex);
      Index = NextIndex;
      MetaOut << FreeListSection;
   }

   MetaOut.flush();
   MetaOut.close();
   FD = -1;

#endif
}

unsigned char const* PageFileImpl::read(size_t Page)
{
   DEBUG_PRECONDITION(!IsOnFreeList(Page))(Page);
   ++PagesRead;

   off_t Offset = off_t(Page) * off_t(PageSize);
   return Alloc->read_file(FD, Offset);
}

size_t PageFileImpl::write(unsigned char const* Buffer)
{
   CHECK(!ReadOnly);

   size_t Page = AllocatePage();

   //   std::cout << "Writing buffer " << (void*) Buffer << " as page " << Page << std::endl;

   ++PagesWritten;

   off_t Offset = off_t(Page) * off_t(PageSize);

#if defined(MULTITHREAD)
   ssize_t Written = ::pwrite(FD, Buffer, PageSize, Offset);
#else
   off_t err = lseek(FD, Offset, SEEK_SET);
   if (err == off_t(-1))
   {
      perror("lseek failed");
      PANIC("fatal")(err)(Offset)(Page)(PageSize);
   }
   ssize_t Written = ::write(FD, Buffer, PageSize);
#endif

   if (Written == -1)
   {
      int Err = errno;
      switch (Err)
      {
         case EBADF  : PANIC("pheap file descriptor is bad.");
         case EFAULT : PANIC("invalid buffer passed to PHeapFileSistem::write().");
         case EFBIG  : PANIC("persistent heap file is too big!");
         case EDQUOT : PANIC("size of persistent heap file exceeds disk quota.");
         case EIO    : PANIC("physical I/O error while writing persistent heap file!");
         case ENOSPC : PANIC("no free space to write persistent heap file.");
         default     : PANIC("cannot write to persistent heap file")(strerror(Err));
      }
   }
   else if (Written < ssize_t(PageSize))
   {
      // partial write
      PANIC("out of disk space while writing persistent heap.")(Written);
   }

   CHECK(Written == ssize_t(PageSize));

   Alloc->deallocate(Buffer);
   return Page;
}

void PageFileImpl::deallocate(size_t Page)
{
   CHECK(Page != 0);  // page 0 is reserved, should never happen
   pthread::mutex::sentry Lock(FreeListMutex);
   DEBUG_PRECONDITION(!IsOnFreeList(Page));
   //   std::cout << "deallcating page " << Page << std::endl;
   FreeList.insert(Page);

   // optimize the free list by removing pages from the free list that are at the end of the file
   while (FreeList.count(NumAllocatedPages-1))
   {
      FreeList.erase(--NumAllocatedPages);
   }
}

size_t PageFileImpl::try_defragment(size_t Page)
{
   size_t FreeListSize;
   {
      pthread::mutex::sentry Lock(FreeListMutex);
      DEBUG_PRECONDITION(!IsOnFreeList(Page));
      FreeListSize = FreeList.size();
   }
   if (Page >= NumAllocatedPages-FreeListSize)
   {
      unsigned char const* Buf = this->read(Page);
      this->deallocate(Page);
      Page = this->write(Buf);
      //      this->get_allocator()->deallocate(Buf);
      // We don't need to deallocate the buffer: write() does it for us
   }
   return Page;
}

size_t PageFileImpl::AllocatePage()
{
   pthread::mutex::sentry Lock(FreeListMutex);
   if (FreeList.empty())
   {
      if (PageCheckpointLimit > 0 && NumAllocatedPages >= PageCheckpointLimit)
      {
         ProcControl::AsyncCheckpoint(ProcControl::ReturnDiskLimitExceeded,
                                      "Page file has exceeded checkpoint limit.");
      }
      return NumAllocatedPages++;
   }
   // else
   size_t Page = *FreeList.begin();
   FreeList.erase(FreeList.begin());
   //   TRACE("Allocated page")(Page);
   CHECK(Page != 0);  // page 0 is reserved, should never happen
   return Page;
}

bool PageFileImpl::IsOnFreeList(size_t Page)
{
   return FreeList.count(Page) == 1;
}

void PageFileImpl::set_checkpoint_limit_kb(unsigned long Size)
{
   PageCheckpointLimit = Size / (PageSize / 1024);
   if (PageCheckpointLimit == 0 && Size > 0) PageCheckpointLimit = 1;
}

unsigned long PageFileImpl::get_checkpoint_limit_kb() const
{
   return PageCheckpointLimit * (PageSize / 1024);
}

void PageFileImpl::Debug()
{
   std::cerr << "PageFileImpl: FileName=" << FileName
             << ", FreeList=" << FreeList
             << '\n';
}

} // namespace PHeapFileSystem
