// -*- C++ -*- $Id$

/*
  2002-07-26: Fixed a bug where FileSystem::flush() would re-write pages if they were in the cache
              but also on disk, resulting in over-use of the page file.
*/

#include "pheapfsv4.h"
#include "rawpagestream.h"
#include "common/inttype.h"
#include <algorithm>
#include <numeric>
#include <errno.h>
#include <string.h>

#if defined(PHEAP_TRACE_DETAILED)
#define TRACE_PHEAP(Msg) TRACE(Msg)
#else
#define TRACE_PHEAP(Msg) DUMMY_TRACE(Msg)
#endif
#if defined(PHEAP_TRACE_DETAILED_EXTRA)
#define TRACE_PHEAP_X(Msg) TRACE(Msg)
#else
#define TRACE_PHEAP_X(Msg) DUMMY_TRACE(Msg)
#endif

namespace PHeapFileSystem
{

char const* BinExtension = ".bin";
char const* MetadataExtension = ".meta";

using inttype::int32;
using inttype::int32_t;
using inttype::uint32;
using inttype::uint32_t;

/*
  Metadata format for FileSystem (XDR): version 1

  Global for each page file:
    byte8 * 4     Magic ("PHFS")
    uint32        Metadata version number
    uint32        Number of page files
    string*N      File names
    uint32        Hook page number
    uint32        PageListSize

  This is immediately followed by data specific to each page file:
    uint32                    file number
    uint32                    Number of allocated pages for this file
    pair(uint32,uint32)*N     For each page, a pair (GlobalPageNumber, LocalPageNumber)

  This is written to a rawpagestream, starting from a page number which is stored in the UserData field
  of each page file.

  The Hook page number is a local page number in the first page file which is the start
  of a rawpagestream which is accessible as metadata_in() or metadata_out().

  Version 2 metadata is identical.  But the changed version number (on 2015/09/16) is
  a hack to allow versioning of objects that don't have their own version number.  The intent
  is that if a file has FileSystem version 1, then it is an old file with a known version number.
  Otherwise if the FileSystem is version 2, then it is a newer file, which is a good
  opportunity to add versioning support to the object!
  */

uint32_t MetadataVersion = 2;

namespace Private
{

//
// PageInfoType
//

PageInfoType::PageInfoType(int InitialLockCount, int InitialReferenceCount, 
			   size_t Page_, unsigned char* Ptr_, FileSystem* FS_, PageFile* PF_, size_t LocalPage_)
   : LockCount(InitialLockCount), ReferenceCount(InitialReferenceCount),
     GlobalPage(Page_), Ptr(Ptr_), BackFS(FS_), PF(PF_), LocalPage(LocalPage_), InPageCache(false)
{
}

void PageInfoType::ResetBuffer(unsigned char* Buf)
{
   DEBUG_PRECONDITION(LockCount.is_zero() && ReferenceCount.is_zero());
   DEBUG_PRECONDITION(Ptr == NULL);
   DEBUG_PRECONDITION(PF == NULL);

   //   DEBUG_TRACE("Resetting PageInfo")(this);

   ++LockCount;
   Ptr = Buf;
   InPageCache = false;
}

void PageInfoType::AddLock()
{
   DEBUG_PRECONDITION(!LockCount.is_zero());
   DEBUG_PRECONDITION(Ptr != NULL);
   ++LockCount;
}

unsigned char* PageInfoType::LoadAddLock()
{
   // 2003-11-16: fixed a multithread race.  Previously, if the LockCount was non-zero then
   // we simply incremented it and returned without ever grabbing the PageMutex (double-checked lock).
   // This is a race because in between the check LockCount.is_zero() and the increment ++LockCount
   // another thread could call SubLock() resulting in the page being removed from memory.
   // This race cannot happen with AddLock() because there it is assumed that the calling
   // thread already has at least one lock.

   // 2007-07-18: There is still a multithread race?  If thread 2 calls SubLock to reduce the lock
   // count to zero, but thread 1 calls LoadAddLock() after thread 2 decrements the counter
   // but before grabbing the mutex, the page might still be in memory.  In this case, we can
   // just do nothing - since we increase the lock count inside the mutex, the double checked lock
   // inside DoLockCountZero will notice.

   pthread::mutex::sentry MyLock(PageMutex);

   TRACE_PHEAP("LoadAddLock")(this)(LocalPage);

   if (LockCount.is_zero())
   {
      TRACE_PHEAP("Calling read()")(this)(LocalPage)(InPageCache);
      BackFS->read(this);
   }
   ++LockCount;
   return Ptr;
}

void PageInfoType::SubLock()
{
   DEBUG_PRECONDITION(!LockCount.is_zero());
   TRACE_PHEAP("SubLock")(this)(LocalPage);
   if (--LockCount == 0) DoLockCountZero();
}

void PageInfoType::AddReference()
{
   ++ReferenceCount;
}

void PageInfoType::SubReference()
{
   if (--ReferenceCount == 0) DoReferenceCountZero();
}

void PageInfoType::DoLockCountZero()
{
   // There is a race condition here.  We allow relaxed decrementing of the LockCount
   // (ie. not with the PageMutex held) in exchange for having to double-check
   // the lock count.  This means it is possible for more than one thread to enter
   // this function when the lock count is zero.  So whatever we do (deallocate, write, whatever),
   // we had better only to it once!

   //   std::cout << "DoLockCountZero PageInfo " << (void*)this << std::endl;
   pthread::mutex::sentry MyLock(PageMutex);

   // re-check the state of LockCount inside the mutex
   if (!LockCount.is_zero()) return;

   TRACE_PHEAP("DoLockCountZero")(this)(LocalPage);
   
   if (BackFS != NULL)
   {
      if (ReferenceCount.is_zero())
      {
	 BackFS->deallocate(this);
      }
      else if (!InPageCache) // if we have entered this function twice, don't write it a second time
      {
	 BackFS->write(this);
      }
   }
}

void PageInfoType::DoReferenceCountZero()
{
   TRACE_PHEAP("DoReferenceCountZero")(this)(LocalPage);
   //   std::cout << "DoReferenceCountZero PageInfo " << (void*)this << std::endl;
   pthread::mutex::sentry MyLock(PageMutex);

   if (ReferenceCount.is_zero() && LockCount.is_zero() && BackFS != NULL)
   {
      BackFS->deallocate(this);
   }
}

} // namespace Private

//
// PageId
//

PageId::PageId() : pInfo(NULL)
{
}

PageId::~PageId()
{
   if (pInfo) pInfo->SubReference();
}

PageId::PageId(Private::PageInfoType* PageInfo) : pInfo(PageInfo)
{
   if (pInfo) pInfo->AddReference();
}

PageId::PageId(PageId const& p) : pInfo(p.pInfo)
{
   if (pInfo) pInfo->AddReference();
}

PageId& PageId::operator=(PageId const& p)
{
   if (p.pInfo) p.pInfo->AddReference();
   if (pInfo) pInfo->SubReference();
   pInfo = p.pInfo;
   return *this;
}

//
// WriteBuffer
//

WriteBuffer::WriteBuffer() : pInfo(NULL), Buf(NULL)
{
}

WriteBuffer::WriteBuffer(Private::PageInfoType* pInfo_, unsigned char* Buf_)
   : pInfo(pInfo_), Buf(Buf_)
{
}

WriteBuffer::~WriteBuffer()
{
   if (pInfo) pInfo->SubLock();
}

WriteBuffer::WriteBuffer(WriteBuffer const& wb)
   : pInfo(wb.pInfo), Buf(wb.Buf)
{
   if (pInfo) pInfo->AddLock();
}

WriteBuffer& WriteBuffer::operator=(WriteBuffer const& wb)
{
   if (wb.pInfo) wb.pInfo->AddLock();
   if (pInfo) pInfo->SubLock();
   pInfo = wb.pInfo;
   return *this;
}

PageId WriteBuffer::get_page() const
{
   return PageId(pInfo);
}

//
// ReadBuffer
//

ReadBuffer::ReadBuffer() : Buf(NULL), pInfo(NULL)
{
}

ReadBuffer::~ReadBuffer()
{
   if (pInfo) pInfo->SubLock();
}

ReadBuffer::ReadBuffer(ReadBuffer const& rb) : Buf(rb.Buf), pInfo(rb.pInfo)
{
   if (pInfo) pInfo->AddLock();
}

ReadBuffer::ReadBuffer(PageId const& Page)
   : pInfo(Page.pInfo)
{
   Buf = pInfo ? pInfo->LoadAddLock() : NULL;
}

ReadBuffer& ReadBuffer::operator=(ReadBuffer const& rb)
{
   if (rb.pInfo) rb.pInfo->AddLock();
   if (pInfo) pInfo->SubLock();
   pInfo = rb.pInfo;
   Buf = rb.Buf;
   return *this;
}

PageId ReadBuffer::get_page() const
{
   return PageId(pInfo);
}

//
// FileSystem
//

// public interface

FileSystem::FileSystem()
   : CurrentPageFile(0),
     PageSize(0),
     HookPageNumber(uint32_t(-1)),
     Alloc(NULL),
     IsReadOnly(false),
     MetaVersion(MetadataVersion)
{
}

FileSystem::~FileSystem()
{
}

void FileSystem::create(std::string const& FilePath, int NumFiles, 
                        size_t DesiredPageSize, size_t PageCacheByteSize,
                        bool Unlink, bool AllowOverwrite)
{
   Alloc = BufferAllocator::GetAllocator(DesiredPageSize);
   PageSize = Alloc->get_page_size();

   std::string::const_iterator PathSep = FilePath.begin();
   std::string::const_iterator PathSepNext = std::find(PathSep, FilePath.end(), '/');
   while (PathSepNext != FilePath.end())
   {
      PathSep = PathSepNext;
      ++PathSep;
      PathSepNext = std::find(PathSep, FilePath.end(), '/');   
   }
   std::string Path = std::string(FilePath.begin(), PathSep);
   if (Path.empty()) Path = ".";  // default path is current directory
   FileName = std::string(PathSep, FilePath.end());

   IsReadOnly = false;
   MaxPageCacheSize = PageCacheByteSize / PageSize + 1;

   for (int i = 0; i < NumFiles; ++i)
   {
      PageFile* PF = new PageFile;
      PF->create(PageSize, Path + '/' + GetFileName(NumFiles, i), Unlink, AllowOverwrite);
      PageFileList.push_back(PF);
      PageFileNames.push_back(GetFileName(NumFiles, i));
      PageFileMetaPages.push_back(std::list<size_t>());
   }
   notify_log(10, pheap::PHeapLog) << "raw file system initialized.  PageSize = " << PageSize << '\n';
}

int FileSystem::ReadPageFileMetadata(std::string const& Path, std::string const& FileName)
{
   notify_log(30, pheap::PHeapLog) << "Reading page file metadata for file " << FileName << '\n';

   // Construct the PageFile
   PageFile* PF = new PageFile;
   uint64 MetaPageNum = PF->open(Path + '/' + FileName, IsReadOnly);

   // check the page size
   size_t CheckPageSize = PF->get_page_size();
   if (PageSize == 0)
   {
      PageSize = CheckPageSize;
   }
   else if (PageSize != CheckPageSize)
   {
      PANIC("Page size mismatch")(FileName);
   }

   if (!Alloc) Alloc = BufferAllocator::GetAllocator(PageSize);

   // read the metadata
   irawpagestream MetaIn(PF, MetaPageNum, PStream::format::XDR);

   // Version number
   char Magic[5];
   Magic[0] = MetaIn.read<char>();
   Magic[1] = MetaIn.read<char>();
   Magic[2] = MetaIn.read<char>();
   Magic[3] = MetaIn.read<char>();
   Magic[4] = 0;

   if (std::string(Magic) != "PHFS")
   {
      PANIC("Invalid page, probably caused by a corrupt file")(FileName);
   }

   uint32_t CheckVersion = MetaIn.read<uint32>();
   MetaVersion = CheckVersion;
   if (CheckVersion < 1 || CheckVersion > 2)
   {
      WARNING("Page file version mismatch, file version is ") << CheckVersion
							      << " expected version is " << MetadataVersion << "\nFilename: " << FileName;
   }

   // Number of page files
   uint32_t CheckNumPageFiles = MetaIn.read<uint32>();
   if (PageFileNames.empty())
   {
      PageFileNames.resize(CheckNumPageFiles);
      PageFileList.resize(CheckNumPageFiles);
   }
   else if (CheckNumPageFiles != PageFileNames.size())
   {
      PANIC("Page file consistency error") << FileName << " lists " 
	    << CheckNumPageFiles << " page files, but we expected " << PageFileNames.size();
   }

   // Page file names
   for (size_t i = 0; i < PageFileNames.size(); ++i)
   {
      std::string CheckFileName = MetaIn.read<std::string>();
      if (PageFileNames[i].empty())
      {
	 PageFileNames[i] = CheckFileName;
      }
      else if (PageFileNames[i] != CheckFileName)
      {
	 PANIC("Page file consistency error") << FileName << " lists a file name of "
	       << CheckFileName << ", but we expected " << PageFileNames[i];
      }
   }

   // HookPageNumber
   uint32_t CheckHookPageNumber = MetaIn.read<uint32>();
   if (HookPageNumber == uint32_t(-1))
   {
      HookPageNumber = CheckHookPageNumber;
   }
   else if (HookPageNumber != CheckHookPageNumber)
   {
      PANIC("Page file consistency error") << FileName << " lists a HookPageNumber of " 
	    << CheckHookPageNumber << " , but we expected " << HookPageNumber;
   }

   // PageListSize
   uint32_t CheckPageListSize = MetaIn.read<uint32>();
   if (PageList.empty())
   {
      PageList.resize(CheckPageListSize);
   }
   else if (PageList.size() != CheckPageListSize)
   {
      PANIC("Page file consistency error") << FileName << " lists a page list size of " 
	    << CheckPageListSize << " , but we expected " << PageList.size();
   }

   // page file number
   uint32_t PageFileNumber = MetaIn.read<uint32>();
   if (PageFileNumber > CheckNumPageFiles)
   {
      PANIC("Page file consistency error") << FileName << " claims to be page file number "
	    << PageFileNumber << ", but only " << CheckNumPageFiles << " files were expected.";
   }

   // Make sure that this PageFileList entry does not yet exist, and then set it
   if (PageFileList[PageFileNumber] != NULL)
   {
      PANIC("Page file consistency error") << FileName << " claims to be file number "
	    << PageFileNumber << ", but this file has already been processed (as file "
	    << PageFileNames[PageFileNumber] << ")";
   }
   PageFileList[PageFileNumber] = PF;

   // If the file names don't match, this is not fatal - but if we have more than 1 file,
   // most likely we will die later trying to read a mis-named file.
   if (PageFileNames[PageFileNumber] != FileName && CheckNumPageFiles > 1)
   {
      WARNING("Page file ") << FileName << " was expected to be named " << PageFileNames[PageFileNumber];
      PageFileNames[PageFileNumber] = FileName;  // assume the supplied filename is authoritative
   }

   // Read the local page information
   uint32_t NumLocalPages = MetaIn.read<uint32>();
   for (uint32_t i = 0; i < NumLocalPages; ++i)
   {
      uint32_t GlobalPageNumber = MetaIn.read<uint32>();
      uint32_t LocalPageNumber = MetaIn.read<uint32>();
      CHECK(PageList[GlobalPageNumber] == NULL);

      PageList[GlobalPageNumber] = new 
	Private::PageInfoType(0, 0, GlobalPageNumber, NULL, this, PF, LocalPageNumber);
   }

   // We no longer free the metadata, so that we don't overwrite it
   // to allow transactions.  Instead we do a deferred free, and store
   // the list of pages used for metadata
   if (PageFileMetaPages.empty())
      PageFileMetaPages.resize(PageFileNames.size());

   PageFileMetaPages[PageFileNumber] = MetaIn.defer_free();

   return PageFileNumber;
}

PageId FileSystem::open(std::string const& FilePath, size_t PageCacheByteSize, bool IsReadOnly_)
{
   notify_log(10, pheap::PHeapLog) << "file system persistent restart initiated\n";
   notify_log(20, pheap::PHeapLog) << "reading raw file system metadata...\n";

   std::string::const_iterator PathSep = FilePath.begin();
   std::string::const_iterator PathSepNext = std::find(PathSep, FilePath.end(), '/');
   while (PathSepNext != FilePath.end())
   {
      PathSep = PathSepNext;
      ++PathSep;
      PathSepNext = std::find(PathSep, FilePath.end(), '/');   
   }
   std::string Path = std::string(FilePath.begin(), PathSep);
   if (*PathSep == '/') ++PathSep;
   if (Path.empty()) Path = ".";  // default path is current directory
   FileName = std::string(PathSep, FilePath.end());

   IsReadOnly = IsReadOnly_;

   int InitialPageFileNumber = ReadPageFileMetadata(Path, FileName);

   for (int FileNum = 0; FileNum < int(PageFileList.size()); ++FileNum)
   {
      if (FileNum == InitialPageFileNumber) continue;         // don't process the initial file twice
      ReadPageFileMetadata(Path, PageFileNames[FileNum]);
   }

   MaxPageCacheSize = PageCacheByteSize / PageSize + 1;

   notify_log(20, pheap::PHeapLog) << "finished reading raw file system metadata:\n";
   notify_log(20, pheap::PHeapLog) << "  PageSize = " << PageSize << '\n';
   notify_log(20, pheap::PHeapLog) << "  NumPageFiles = " << PageFileList.size() << '\n';
   notify_log(20, pheap::PHeapLog) << "  PageListSize = " << PageList.size() << '\n';
   notify_log(20, pheap::PHeapLog) << "  FreeListSize = " << FreeList.size() << '\n';

   return HookPageNumber == uint32_t(-1) ? PageId() : this->get_page(HookPageNumber);
}

void FileSystem::defragment()
{
   notify_log(20, pheap::PHeapLog) << "Defragmenting page files...\n";
   pthread::mutex::sentry MyLock(PageCacheMutex);

   for (std::size_t i = 0; i < PageFileList.size(); ++i)
   {
      while (PageFileList[i]->num_free_pages() > 0)
      {
         // find a page to reallocate
         for (size_t page = 0; page < PageList.size(); ++page)
         {
            if (PageList[page] && PageList[page]->PF == PageFileList[i])
               PageList[page]->LocalPage = PageList[page]->PF->try_defragment(PageList[page]->LocalPage);
         }
      }
   }
}

void FileSystem::flush()
{
   notify_log(20, pheap::PHeapLog) << "Flushing page files...\n";
   pthread::mutex::sentry MyLock(PageCacheMutex);

   // firstly, clean out the page cache
   while (!PageCache.empty())
   {
      Private::PageInfoType* PageInfo = PageCache.back();
      PageCache.pop_back();

      //      DEBUG_TRACE("Removed PageInfo from PageCache")((void*) PageInfo)(PageCache.size());

      CHECK(PageInfo->InPageCache);
      CHECK(PageInfo->Ptr != NULL);

      if (PageInfo->PF == NULL)
      {
	 PageInfo->PF = AllocatePageFile();
	 PageInfo->LocalPage = PageInfo->PF->write(PageInfo->Ptr);
      }
      else
      {
	 DeallocatePageBuffer(PageInfo->Ptr);
      }
      PageInfo->InPageCache = false;
      PageInfo->Ptr = NULL;
   }

   // Now flush all remaining pages to disk if they are in memory but not yet on disk.
   // This implies that the lock count is non-zero.  We should be able to avoid writing
   // pages that have a non-zero lock count but zero reference count, but this is not
   // yet done.
   for (size_t i = 0; i < PageList.size(); ++i)
   {
      if (PageList[i] && PageList[i]->Ptr != NULL && PageList[i]->PF == NULL)
      {
	 CHECK(!PageList[i]->LockCount.is_zero());
	 PageList[i]->PF = AllocatePageFile();
	 PageList[i]->LocalPage = PageList[i]->PF->write(PageList[i]->Ptr);
      }
   }
   notify_log(20, pheap::PHeapLog) << "Flushing page files completed.\n";
}	 

void FileSystem::persistent_shutdown(PageId UserPage)
{
   // Now we can delete the old metadata pages that we deferred deallocating
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      while (!PageFileMetaPages[i].empty())
      {
         PageFileList[i]->deallocate(PageFileMetaPages[i].front());
         PageFileMetaPages[i].pop_front();
      }
   }

   // Defragment pages and flush all in-memory pages to disk
   this->defragment();
   this->flush();

   notify_log(20, pheap::PHeapLog) << "writing pheapfsv4 metadata...\n";

   // Write our own metadata
   if (UserPage.pInfo) UserPage.pInfo->AddReference();  // Add a ref count for the UserPage
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      //      std::cout << "Writing page metadata for file " << i << std::endl;
      this->WritePageFileMetadata(i, UserPage.pInfo ? UserPage.pInfo->GlobalPage : 0);
   }


   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      delete PageFileList[i];
   }
   PageFileList.clear();
   PageFileNames.clear();

   notify_log(10, pheap::PHeapLog) << "file system shutdown completed.\n";
   notify_log(10, pheap::PHeapLog) << "Total pages allocated: " << num_allocated_pages() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total file pages allocated: " << num_file_allocated_pages() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total pages free: " << num_free_pages() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total pages read: " << num_pages_read() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total pages written: " << num_pages_written() << '\n';
}

void FileSystem::WritePageFileMetadata(int FileNumber, uint32_t MetaHook)
{
   PageFile* PF = PageFileList[FileNumber]; // this is the page file that we are concerned with
   // allocate a rawpagestream to write the metadata
   orawpagestream out(PF, PStream::format::XDR);

   // this data is the same for every file
   char const* Magic = "PHFS";
   out.write<char>(Magic[0]);
   out.write<char>(Magic[1]);
   out.write<char>(Magic[2]);
   out.write<char>(Magic[3]);
   out.write<uint32>(MetadataVersion);
   out.write<uint32>(PageFileList.size());
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      out << PageFileNames[i];
   }
   out.write<uint32>(MetaHook);
   out.write<uint32>(PageList.size());

   // this data is specific to this file number
   out.write<uint32>(FileNumber);
   // count how many pages are allocated to this file
   size_t NumLocalPages = 0;
   for (size_t i = 0; i < PageList.size(); ++i)
   {
      if (PageList[i] && PageList[i]->PF == PF) ++NumLocalPages;
   }

   out.write<uint32>(NumLocalPages);
   // write the pages
   for (size_t i = 0; i < PageList.size(); ++i)
   {
      if (PageList[i] && PageList[i]->PF == PF)
      {
	 out.write<uint32>(PageList[i]->GlobalPage);
	 out.write<uint32>(PageList[i]->LocalPage);
      }
   }
   uint32_t LocalHook = out.commit();
   PF->persistent_shutdown(LocalHook);
}

WriteBuffer FileSystem::allocate()
{
   pthread::mutex::sentry MyLock(PageCacheMutex);

   unsigned char* Buf = AllocatePageBuffer();
   Private::PageInfoType* Page;
   if (FreeList.empty())
   {
      Page = new Private::PageInfoType(1, 0, PageList.size(), Buf, this, NULL, 0);
      //      std::cout << "Allocated new PageInfo " << Page->Page << ' ' 
      //		<< PageList.size() << ' ' << FreeList.size() << ' ' << num_free_pages() << ' '
      //		<< num_file_allocated_pages() << std::endl;
      PageList.push_back(Page);
   }
   else
   {
      std::size_t FreePage = *FreeList.begin();
      FreeList.erase(FreeList.begin());
      Page = PageList[FreePage];
      Page->ResetBuffer(Buf);
      //      std::cout << "Allocated from free list PageInfo " << (void*) Page << ' ' 
      //		<< PageList.size() << ' ' << FreeList.size() << ' ' << num_free_pages() << ' '
      //		<< num_file_allocated_pages() << std::endl;
   }
   return WriteBuffer(Page, Buf);
}

PageId FileSystem::get_page(size_t Page) const
{
   // getting the PageCacheMutex here is perhaps a little over cautious; this function should
   // only be called when restarting a saved state, via the PageId extractor.  It is unlikely
   // (and possibly a bug) if we allocate a page or otherwise do something that may invalidate
   // the PageList.  However it is technically correct.
   Private::PageInfoType* PageInfo;
   {
      pthread::mutex::sentry MyLock(PageCacheMutex);
      PageInfo = PageList[Page];
   }
      
   // Make sure the page exists.  
   // Technically, we should grab the PageMutex to do the precondition check...
   DEBUG_PRECONDITION(PageInfo->Ptr != NULL || PageInfo->PF != NULL);
   return PageId(PageInfo); 
}



// private interface for PageInfoType.  The caller holds PageInfo->PageMutex for these functions.

PageFile* FileSystem::AllocatePageFile()
{
   PageFile* Temp = PageFileList[CurrentPageFile];
   CurrentPageFile = (CurrentPageFile+1) % PageFileList.size();
   //   std::cout << "allocated PageFile " << Temp->Page << std::endl;
   return Temp;
}

void FileSystem::read(Private::PageInfoType* PageInfo)
{
   //   std::cout << "reading PageInfo " << PageInfo->Page << std::endl;
   // get the PageCacheMutex and see if its in the page cache
   {
      if (PageInfo->InPageCache)
      {
         pthread::mutex::sentry MyLock(PageCacheMutex);
	 //      	 std::cout << "reading PageInfo " << PageInfo->Page << " is in cache." << std::endl;
	 PageCache.erase(PageInfo->PageCacheLoc);
	 PageInfo->InPageCache = false;

         //	 DEBUG_TRACE("Removed PageInfo from PageCache")((void*) PageInfo)(PageCache.size());

	 return;
      }
   }

   if (PageInfo->Ptr != NULL)
      return;

   //   std::cout << "reading PageInfo " << PageInfo->Page << " from disk." << std::endl;
   // else not in the page cache.  Read it from disk.
   CHECK(PageInfo->PF != NULL);
   CHECK(PageInfo->Ptr == NULL);
   // the const_cast here is slightly dangerous in that it is possible that we don't have
   // write permission on the buffer.  No issue here though, since we never re-write buffers anyway.
   PageInfo->Ptr = const_cast<unsigned char*>(PageInfo->PF->read(PageInfo->LocalPage));
}

void FileSystem::write(Private::PageInfoType* PageInfo)
{
   pthread::mutex::sentry MyLock(PageCacheMutex);

   TRACE_PHEAP("Writing page")(PageInfo->LocalPage);

   DEBUG_PRECONDITION(PageInfo->LockCount.is_zero() && !PageInfo->ReferenceCount.is_zero());
   DEBUG_PRECONDITION(!PageInfo->InPageCache);
   DEBUG_PRECONDITION(PageInfo->Ptr != NULL);

   PageCache.push_front(PageInfo);
   PageInfo->InPageCache = true;
   PageInfo->PageCacheLoc = PageCache.begin();

   //   DEBUG_TRACE("Added PageInfo to PageCache")((void*) PageInfo)(PageCache.size());

   // if the page cache is too big, flush some pages to disk.
   while (!is_read_only() && PageCache.size() > MaxPageCacheSize)
   {
      // pop a page off the cache
      PageInfo = PageCache.back();
      PageCache.pop_back();

      //      DEBUG_TRACE("Removed PageInfo from PageCache")((void*) PageInfo)(PageCache.size());

      CHECK(PageInfo->InPageCache);
      CHECK(PageInfo->Ptr != NULL);

      //      std::cout << "removing PageInfo " << (void*)PageInfo << " from cache." << std::endl;
      // Only write the page if its not already on disk
      if (PageInfo->PF == NULL)
      {
	 // 	 std::cout << "flushing PageInfo " << PageInfo->Page << " to disk." << std::endl;
	 PageInfo->PF = AllocatePageFile();
	 PageInfo->LocalPage = PageInfo->PF->write(PageInfo->Ptr);
	 //std::cout << "local page is " << PageInfo->LocalPage << std::endl;
      }
      else
      {
	 // The page is already on disk, just drop the in-memory version
	 DeallocatePageBuffer(PageInfo->Ptr);
      }
      PageInfo->InPageCache = false;
      PageInfo->Ptr = NULL;
   }
}

void FileSystem::deallocate(Private::PageInfoType* PageInfo)
{
   TRACE_PHEAP("deallocate")(PageInfo);
   pthread::mutex::sentry MyLock(PageCacheMutex);
   DEBUG_PRECONDITION(PageInfo->LockCount.is_zero() && PageInfo->ReferenceCount.is_zero());
   if (PageInfo->InPageCache)
   {
      PageCache.erase(PageInfo->PageCacheLoc);
      PageInfo->InPageCache = false;

      //      DEBUG_TRACE("Removed PageInfo from PageCache")((void*) PageInfo)(PageCache.size());
   }
   if (PageInfo->Ptr != NULL)
   {
      DeallocatePageBuffer(PageInfo->Ptr);
      PageInfo->Ptr = NULL;
   }
   // If the page is on disk then deallocate it.
   // If the PageFileList is empty, then we have already shutdown, which is probably a bug
   // but we can shutdown successfully
   if (PageInfo->PF != NULL) // && !PageFileList.empty())
   {
      if (PageFileList.empty())
      {
	 notify_log(0, pheap::PHeapLog) 
	   << "probable bug: in FileSystem::deallocate(Private::PageInfoType* PageInfo): "
	   << "deallocating a page but the PageFile no longer exists!\n";
	 return;
      }
      PageInfo->PF->deallocate(PageInfo->LocalPage);
      //      std::cout << "Local page is " << PageInfo->LocalPage << std::endl;
      PageInfo->PF = NULL;
   }
   FreeList.insert(PageInfo->GlobalPage);
}

std::string FileSystem::GetFileName(int NumFiles, int FileNumber) const
{
   std::ostringstream FName;
   FName << FileName;
   if (NumFiles > 1) FName << '.' << (FileNumber+1);
   return FName.str();
}

void FileSystem::shutdown(bool Remove)
{
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      PageFileList[i]->shutdown(Remove);
   }

   PageCache.clear();

   // walk the page list and deallocate the pages.  We don't destroy the PageInfo*'s themselves,
   // since there may be pointers to them elsewhere (as PageId's, for example) that
   // will want to destroy themselves cleanly.  
   for (size_t i = 0; i < PageList.size(); ++i)
   {
      if (PageList[i])
      {
	 PageList[i]->BackFS = NULL;
	 PageList[i]->PF = NULL;
	 PageList[i]->InPageCache = false;
	 if (PageList[i]->Ptr != NULL)
	 {
	    DeallocatePageBuffer(PageList[i]->Ptr);
	    PageList[i]->Ptr = NULL;
	 }
      }
   }

   notify_log(10, pheap::PHeapLog) << "file system shutdown completed.\n";
   notify_log(10, pheap::PHeapLog) << "Total pages allocated: " << num_allocated_pages() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total file pages allocated: " << num_file_allocated_pages() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total pages free: " << num_free_pages() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total pages read: " << num_pages_read() << '\n';
   notify_log(10, pheap::PHeapLog) << "Total pages written: " << num_pages_written() << '\n';
}

int FileSystem::GetPageFileNumber(PageFile* PF) const
{
   int Loc = std::find(PageFileList.begin(), PageFileList.end(), PF) - PageFileList.begin();
   return Loc == int(PageFileList.size()) ? -1 : Loc;
}

size_t FileSystem::num_allocated_pages() const
{
   return PageList.size();
}

int FileSystem::num_pages_written() const
{
   int Total = 0;
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      Total += PageFileList[i]->num_pages_written();
   }
   return Total;
}

int FileSystem::num_file_allocated_pages() const
{
   int Total = 0;
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      Total += PageFileList[i]->num_allocated_pages();
   }
   return Total;
}

int FileSystem::num_free_pages() const
{
   int Total = 0;
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      Total += PageFileList[i]->num_free_pages();
   }
   return Total;
}

int FileSystem::num_pages_read() const
{
   int Total = 0;
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      Total += PageFileList[i]->num_pages_read();
   }
   return Total;
}

void FileSystem::set_checkpoint_limit_kb(unsigned long Size)
{
   notify_log(20, pheap::PHeapLog) << "Setting checkpoint limit to " << Size << " KB\n";
   unsigned long SizePerFile = Size / PageFileList.size();
   if (SizePerFile == 0 && Size > 0) SizePerFile = 1;  // make sure small values of Size don't get rounded to zero
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      PageFileList[i]->set_checkpoint_limit_kb(SizePerFile);
   }
   notify_log(20, pheap::PHeapLog) << "checkpoint limit is " << this->get_checkpoint_limit_kb() << " KB\n";
}

unsigned long FileSystem::get_checkpoint_limit_kb()
{
   unsigned long Total = 0;
   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      Total += PageFileList[i]->get_checkpoint_limit_kb();
   }
   return Total;
}

void FileSystem::write_page_id(PStream::opstream& out, PageId const& pid) const
{
   if (pid.pInfo)
   {
      DEBUG_CHECK(pid.pInfo->GetFS() == this);
      out.write<uint32>(pid.pInfo->GlobalPage);
      pid.pInfo->AddReference();
   }
   else
   {
      out.write<uint32>(uint32_t(-1));
   }
}

PageId FileSystem::read_page_id(PStream::ipstream& in) const
{
   uint32_t Page = in.read<uint32>();
   return Page == uint32_t(-1) ? PageId() : this->get_page(Page);
}

std::ostream& operator<<(std::ostream& out, std::list<std::size_t> const& l)
{
   std::copy(l.begin(), l.end(), std::ostream_iterator<std::size_t>(out, ", "));
   return out;
}

void FileSystem::Debug()
{
   std::cerr << "FileSystem debug report\n";
   //	     << "FreeList: " << FreeList;

   for (size_t i = 0; i < PageFileList.size(); ++i)
   {
      PageFileList[i]->Debug();
   }

   std::cerr << '\n';
}

} // namespace PHeapFileSystem
