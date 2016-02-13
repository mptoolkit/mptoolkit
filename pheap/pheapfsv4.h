/* -*- C++ -*- $Id$

  Version 4 of the persistent heap file system.

  Creaded 2002-07-11 Ian McCulloch

  This version uses a write-once, read-many model.  WriteBuffer's are allocated using
  the FileSystem::allocate() function.  A WriteBuffer encapsulates a writable page in memory.
  Pages are read back via the FileSystem::read(PageId) function.  This returns a buffer
  that is read-only.

  The PageId type is a reference-counted handle to a page.
  WriteBuffer and ReadBuffer maintain a per-page lock count.  When the lock count reaches zero,
  the page is flushed to disk (or possibly put on the cache), and the memory buffer is returned
  to the buffer pool. If both the lock and reference counts go to zero, the page is deallocated.
  It is not necessary to explicitly deallocate pages.

  It would be possible to make a re-writable model.  this would require being able to construct
  a WriteBuffer from a PageId.  This function would have to ensure that the page is in memory,
  and 'remove' it from the disk file, and set pInfo->PF = NULL.

  A FileSystem holds pointers to one or more PageFile's.  Pages are allocated to the PageFile's
  on a round-robin basis.

  Still to be done:
  Possibly make the 'get next page file' function in FileSystem class a virtual function.
  This would allow some cool things, like shadow file-systems etc.
  Some of the PageInfoType & PageId & ReadBuffer & WriteBuffer functions should be inlined.
  There is no locking in the init/shutdown functions.  It might be hard to do this correctly.
  eg, we probably should get the mutex on every _page_, as well as the PageCacheMutex.  
  But that itself isn't really possible (ie, we need the PageCacheMutex to
  access PageList, but we cannot get a lock on a PageMutex with PageList alreeady locked without
  a possible deadlock).  Ultimately, without serializing everything with a single lock,
  the only way to properly shutdown the FS is for the overall caller to be serialized.
*/

#if !defined(PHEAPFSV4_H_VEUIOH8RTY89YG78Y7834R8934Y89EYRE)
#define PHEAPFSV4_H_VEUIOH8RTY89YG78Y7834R8934Y89EYRE

#include "common/atomic.h"
#include "common/mutex.h"
#include "common/inttype.h"
#include "pagefile.h"
#include "pstream/pstreamfwd.h"
#include <vector>
#include <list>

namespace pheap
{
extern MessageLogger::Logger PHeapLog;

// if the ExpectedPageFileVersion is set, then
// abort if the file has the wrong version.
// The default expexted page file version is 2
int ExpectedPageFileVersion();

} // namespace pheap

namespace PHeapFileSystem
{

class FileSystem;
class PageId;
 std::ostream& operator<<(std::ostream& out, PageId const& pid);

namespace Private   // Stuff for internal use
{

class PageInfoType
{
   public:
      PageInfoType(int InitialLockCount, int InitialReferenceCount, 
		   size_t Page_, unsigned char* Ptr_, FileSystem* FS_, PageFile* PF_, size_t LocalPage_);

      void ResetBuffer(unsigned char* Buf); // when reference & lock counts are zero, this resets the page
      // info to a new page, with lock count = 1, reference count = 0, Ptr = Buf.

      unsigned char* LoadAddLock();  // same as AddLock(), but allows for the possibility 
                            // that the lock count was previously zero.  Returns a pointer to the buffer.
      void AddLock();
      void SubLock();

      void AddReference();
      void SubReference();

      // returns the file system associated with the page.  The FS doesn't change once the
      // PageInfoType is constructed, so no locking is required for this function.
      FileSystem* GetFS() const { return BackFS; }

      // writes the Page number to the stream.  This function is not, and cannot, be threadsafe.
      // Any modificiations to the underlying page or FileSystem will render the information
      // written by this function incorrect.
      // The opstream cannot be just any opstream, it must be a metadata stream.
      // The corresponding extractor does not have the same threading issues, and is 
      // part of PageId instead.
      void WriteStream(PStream::opstream& out);

   private:
      PageInfoType(PageInfoType const&); // not implemented
      PageInfoType& operator=(PageInfoType const&); // not implemented

      void DoLockCountZero();  // called by SubLock if the lock count hits zeron
      void DoReferenceCountZero(); // called by SubReference if the reference count hits zero

      typedef std::list<PageInfoType*> PageListType;
      typedef PageListType::iterator PageListIterType;

      AtomicRefCount LockCount;
      AtomicRefCount ReferenceCount;

      size_t GlobalPage;             // Page number tracked by the FileSystem
      unsigned char* Ptr;            // pointer to the page, if it is in memory
      FileSystem* BackFS;            // Backpointer to the owning file system
      PageFile* PF;                  // this is non-NULL if the page is on-disk, NULL otherwise.
      size_t LocalPage;              // Used privately by *PF

      bool InPageCache;              // true if *BackFS has this page in cache, false otherwise.
      PageListIterType PageCacheLoc; // if InPageCache is true, this is a backpointer to the cache location

      pthread::mutex PageMutex;      // This mutex protects all the data (except the atomic counters)

   friend class PHeapFileSystem::FileSystem;
   friend class PHeapFileSystem::PageId;

   friend std::ostream& PHeapFileSystem::operator<<(std::ostream& out, PageId const& pid);

   friend std::ostream& operator<<(std::ostream& out, PageInfoType const& p)
   {
      out << "LocalPage: " << p.LocalPage;
      return out;
   }

};

} // namespace Private

class ReadBuffer;
class WriteBuffer;

class PageId
{
   public:
      PageId();
      ~PageId();
      PageId(PageId const& Page);
      PageId& operator=(PageId const& Page);

      // debug functions to bump up the reference count - do not use!
      void bump_count() { pInfo->AddReference(); }
      void sub_count() { pInfo->SubReference(); }

      FileSystem* get_fs() const { return pInfo->GetFS(); }

      bool is_null() const { return pInfo == NULL; }

   private:
      explicit PageId(Private::PageInfoType* PageInfo);

      Private::PageInfoType* pInfo;

   friend class ReadBuffer;
   friend class WriteBuffer;
   friend class FileSystem;

   friend bool operator==(PageId const& p1, PageId const& p2);
   friend bool operator!=(PageId const& p1, PageId const& p2);

   friend std::ostream& operator<<(std::ostream& out, PageId const& pid)
   {
      out << "(GlobalPage: " << (pid.pInfo ? int32_t(pid.pInfo->GlobalPage) : -1);
      if (pid.pInfo) out << ", " << *pid.pInfo;
      out << ')';
      return out;
   }
};

inline bool operator==(PageId const& p1, PageId const& p2)
{
   return p1.pInfo == p2.pInfo;
}

inline bool operator!=(PageId const& p1, PageId const& p2)
{
   return p1.pInfo != p2.pInfo;
}

// WriteBuffer is a reference counted wrapper around a writable buffer
class WriteBuffer
{
   public:
      WriteBuffer();
      ~WriteBuffer();
      WriteBuffer(WriteBuffer const& wb);
      WriteBuffer& operator=(WriteBuffer const& wb);

      unsigned char* buffer() const { return Buf; }

      PageId get_page() const;

      FileSystem* get_fs() const { return pInfo->GetFS(); }

   private:
      WriteBuffer(Private::PageInfoType* pInfo_, unsigned char* Buf_);
      
      Private::PageInfoType* pInfo;
      unsigned char* Buf;

   friend class FileSystem;
};

class ReadBuffer
{
   public:
      ReadBuffer();
      ~ReadBuffer();
      ReadBuffer(ReadBuffer const& rb);
      ReadBuffer& operator=(ReadBuffer const& rb);

      explicit ReadBuffer(PageId const& Page);

      PageId get_page() const;

      unsigned char const* buffer() const { return Buf; }

      FileSystem* get_fs() const { return pInfo->GetFS(); }

   private:
      unsigned char const* Buf;
      Private::PageInfoType* pInfo;
};

class FileSystem
{
   public:
      FileSystem();
      virtual ~FileSystem();

      void create(std::string const& FileName, int NumFiles, 
                  size_t PageSize, size_t PageCacheByteSize, bool Unlink = false, bool AllowOverwrite = true);

      // close the associated page files.  If Remove is true then delete the files.
      void shutdown(bool Remove = false);

      // if the FileSystem was created (rather than opened from an existing file),
      // then assume that the files are now unusable and delete them.
      void cleanup();

      void persistent_shutdown(PageId UserPage);

      int version() { return MetaVersion; }

      PageId open(std::string const& FileName, size_t PageCacheByteSize, bool ReadOnly = false);

      // indicates that an asynchronous checkpoint should be triggered if the size of allocated pages (in KB)
      // exceeds the given limit.  0 means no limit.
      void set_checkpoint_limit_kb(unsigned long Size);

      // returns the current checkpoint limit
      unsigned long get_checkpoint_limit_kb();

      bool is_read_only() const { return IsReadOnly; }
   
      WriteBuffer allocate();

      size_t get_page_size() const { return PageSize; }

      size_t get_page_cache_size() const { return MaxPageCacheSize; }

      // some statistics
      size_t num_allocated_pages() const;
      int num_file_allocated_pages() const;
      int num_free_pages() const;
      int num_pages_written() const;
      int num_pages_read() const;

      void write_page_id(PStream::opstream& out, PageId const& pid) const;
      PageId read_page_id(PStream::ipstream& in) const;

      void defragment();

      void Debug();

   private:
      std::string GetFileName(int NumFiles, int FileNumber) const;

      unsigned char* AllocatePageBuffer() const { return Alloc->allocate(); }
      void DeallocatePageBuffer(unsigned char const* Buf) const { Alloc->deallocate(Buf); }

      // gets a PageId for the given page.  This is used by the PageId stream extractor.
      // The given page must be allocated (ie not on the free list).
      PageId get_page(size_t Page) const;

      // this interface is used by PageInfoType.  The caller must hold the PageMutex
      // before calling these functions.

      // Reads a page into memory.
      // Precondition: PageInfo->Ptr == NULL.
      void read(Private::PageInfoType* PageInfo);

      // deallocates the given page.
      void deallocate(Private::PageInfoType* PageInfo);

      // writes the page to disk (or puts it in the cache).
      // precondition: PageInfo->Ptr != NULL.
      void write(Private::PageInfoType* PageInfo);

      // This interface is used internally by FileSystem.  The caller must hold the PageCacheMutex,
      // as well as the PageMutex for the appropriate page.

      // Returns a PageFile* that is the next target for a written page.
      // Currently, these are allocated round-robin among the possible PageFiles.
      // More sophisticated algorithms would be possible.
      PageFile* AllocatePageFile();

      // Flushes all in-memory pages to disk.  This is not a public function as it
      // is a race condition to write a page to disk if the lock count is > 0.
      // It is only 'safe' to do this immediately before shutting down, when we 'know'
      // that no pages currently in memory will be written to.  
      void flush();

      // utility function to convert a PageFile* into an integer index into the
      // PageFileList.  Returns -1 if PF is not in the PageFileList.
      int GetPageFileNumber(PageFile* PF) const;

      int ReadPageFileMetadata(std::string const& Path, std::string const& FileName);
      void WritePageFileMetadata(int FileNumber, inttype::uint32_t MetaHook);

      mutable pthread::mutex PageCacheMutex;

      size_t MaxPageCacheSize; // TODO: This is not set anywhere

      int CurrentPageFile;  // an index into the PageFileList.
      std::vector<PageFile*> PageFileList;
      std::vector<std::string> PageFileNames;
      std::vector<std::list<size_t> > PageFileMetaPages;

      size_t PageSize;  // size in bytes of each page

      inttype::uint32_t HookPageNumber;  // the page number to be returned by open()

      std::string FileName;

      std::vector<Private::PageInfoType*> PageList;

      int PagePageCacheSize;
      std::list<Private::PageInfoType*> PageCache;

      typedef std::set<size_t> FreeListType;
      FreeListType FreeList;

      BufferAllocator* Alloc;

      bool CreatedNew;  // set to true if the files were created,
                        // so that on an abnormal shutdown we can safely delete them
      bool IsReadOnly;

      int MetaVersion;

   friend class Private::PageInfoType;

   // friendship of PageInfoType isn't transitive to the extractor, it needs to be explicit
   friend PStream::ipstream& operator>>(PStream::ipstream& in, PageId& pid);
};

} // namespace PHeapFileSystem

#endif
