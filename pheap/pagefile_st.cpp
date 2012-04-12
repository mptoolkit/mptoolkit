// -*- C++ -*- $Id$

/*
  pagefile_st.cpp

  Single-threaded implementation of the page file interface.
  Although this is single-threaded, it is thread safe and reenterable.  maybe.

  Created in antiquity, Ian McCulloch
*/

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
  the format of the per-file metadata is:

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
*/

namespace PHeapFileSystem
{

using PStream::ipfilestream;
using PStream::opfilestream;

using inttype::uint32;
using inttype::uint64;

// a 64-bit 'random' number that identifies binary page files.
uint64 PageFileMagic = 0x1fe03dc25ba47986LL;
uint32 PageFileMetadataVersion = 1;

class PageFileImpl
{
   public:
      PageFileImpl();

   void create(size_t PageSize_, std::string const& FileName_, bool Unlink = false);

      uint64 open(std::string const& FileName_, bool ReadOnly);

      void shutdown(bool Remove = false);

      void persistent_shutdown(uint64 UserData);

      size_t write(unsigned char const* Buffer);

      unsigned char const* read(size_t Page);

      void deallocate(size_t Page);

      size_t try_defragment(size_t Page);

      BufferAllocator* get_allocator() const { return Alloc; }

      std::string const& name() const { return FileName; }

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

      pthread::mutex FreeListMutex;
      std::set<size_t> FreeList;

      bool ReadOnly;

      int PagesRead, PagesWritten;

      unsigned long PageCheckpointLimit;
};

PageFileImpl::PageFileImpl()
  : Alloc(NULL), FD(-1), PageSize(0), NumAllocatedPages(0), ReadOnly(false),
     PagesRead(0), PagesWritten(0), 
     PageCheckpointLimit(0)
{
}

void PageFileImpl::create(size_t PageSize_, std::string const& FileName_, bool Unlink)
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

   notify_log(20, pheap::PHeapLog) << "creating file " << FileName << '\n';
   FD = ::open(FileName.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0600);
   if (FD == -1)
   {
      PANIC("Error creating page file!")(FileName)(strerror(errno));
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

   FD = ::open(FileName.c_str(), OpenFlags, 0666);
   if (FD == -1)
   {
      PANIC("Error opening page file!")(FileName)(strerror(errno));
   }

   // create a ipfilestream to read the metadata
   ipfilestream MetaIn(PStream::format::XDR);
   MetaIn.set_fd(FD);
   uint64 Magic = MetaIn.read<uint64>();

   if (Magic != PageFileMagic)
   {
     PANIC(FileName_ + " does not appear to be a page file!")(Magic)(PageFileMagic);
     // << " Magic = " 
     //	   << std::hex << Magic << ", expected = " << std::hex << PageFileMagic;
   }

   uint32 Version = MetaIn.read<uint32>();
   notify_log(40, pheap::PHeapLog) << "PageFile " << FileName_ << " version number is " << Version << '\n';

   if (Version != 1)
   {
      PANIC("PageFile version mismatch")(FileName_)(Version) << "expected version 1";
   }

   PageSize = MetaIn.read<uint32>();
   NumAllocatedPages = MetaIn.read<uint32>();
   uint64 UserData = MetaIn.read<uint64>();

   // Now we know the page size, get an allocator
   Alloc = BufferAllocator::GetAllocator(PageSize);

   // retrieve the free list.  This is a vector of uint32, we need
   // to convert from the default size_t
   off_t Offset = off_t(NumAllocatedPages) * off_t(PageSize);
   MetaIn.lseek(Offset, SEEK_SET);

   std::vector<inttype::uint32> TempFreeList;
   MetaIn >> TempFreeList;

   FreeList = std::set<size_t>(TempFreeList.begin(), TempFreeList.end());

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
}

unsigned char const* PageFileImpl::read(size_t Page)
{
   DEBUG_PRECONDITION(!IsOnFreeList(Page))(Page);
   ++PagesRead;

   off_t Offset = off_t(Page) * off_t(PageSize);

   return Alloc->read_file(FD, Offset);
}

#if 0 // old implementation, before BufferAllocator::read_file() existed

   unsigned char* Buffer = Alloc->allocate();
   //   std::cout << "Reading page " << Page << " into buffer " << (void*) Buffer << std::endl;
   ssize_t Read = ::pread(FD, Buffer, PageSize, Offset);

   if (Read == -1)
   {
      int Err = errno;
      switch (Err)
      {
         case EBADF  : PANIC("pheap file descriptor is bad.");
         default     : PANIC("cannot read from persistent heap file")(Err);
      }
   }
   else if (Read < PageSize)
   {
      // partial read
      PANIC("Did not read a complete page - is the file corrupt?");
   }

   return Buffer;
}
#endif

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
      this->get_allocator()->deallocate(Buf);
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
