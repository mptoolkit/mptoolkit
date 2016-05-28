// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheap.h
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
  Yet another version of the persistent heap.

  Created 2002-12-04 Ian McCulloch
*/

#if !defined(PHEAP_H_HUIHFE89Y54T98YA89TATP5ATO)
#define PHEAP_H_HUIHFE89Y54T98YA89TATP5ATO

#include "common/trace.h"
#include "common/atomic.h"
#include "common/poolallocator.h"
#include "common/mutex.h"
#include "common/hash_map.h"
#include "common/niftycounter.h"
#include "pstream/pstream.h"
#include "pheaperror.h"
#include "pheapallocator.h"
#include "pheapfsblock.h"
#include "pheapstream.h"  // TODO: we shouldn't really need this, only for ::Write() function
#include <list>
#include <typeinfo>

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

namespace pheap
{

using PStream::id_type;

namespace Private
{

void InitializePHeap();

namespace
{
NiftyCounter::nifty_counter<InitializePHeap> InitPHeap;
} // namespace

using namespace PHeapFileSystem;

class PointerToObject
{
   public:
      virtual ~PointerToObject() {}
      virtual PointerToObject* clone() const = 0;

      virtual void WriteToStream(PStream::opstream& out) const = 0;

      virtual std::type_info const& GetTypeid() const = 0;

      template <typename T>
      static PointerToObject* Create(T* Ptr);

      template <typename T>
      T* GetPointer() const;

   protected:
      PointerToObject() {}

   private:
      PointerToObject(const PointerToObject&);             // not defined
      PointerToObject& operator=(const PointerToObject&);  // not defined
};

} // namespace Private

//
// Replacement of the old PHeapObject.  It maintains
// its own reference count, and will self-destruct when the reference count hits zero.
// No public constructors, use the PHeapObject::Create() functions.
//

class PHeapObject
{
   public:
      // constructs a PHeapObject from the given object.
      // The inital LockCount is 1.
      template <typename T>
      static PHeapObject* Create(T* Obj);

      // constructs a PHeapObject from the given object,
      // using the specified ObjectID.
      // For use when importing an object that already has an ID.
      template <typename T>
      static PHeapObject* Create(T* Obj, id_type ObjectID_);

      static PHeapObject* Create(Private::PointerToObject* Obj, id_type ObjectID_);

      // creates a copy of a descriptor (ie from another filesystem).
      static PHeapObject* Create(PHeapFileSystem::Descriptor const& Desc, 
				 id_type ID, 
				 int InitialReferenceCount);

      // creates a PHeapObject with no Object or Descriptor,
      // and with the mutex locked.  This adds the object to the heap,
      // but it cannot be accessed until FinalizeCreate() is called
      // to unlock the mutex.  The initial reference count is 1.
      // This function *ASSUMES THAT GlobalHeapMutex IS LOCKED BY THE CALLER*
      static PHeapObject* CreateLocked(id_type ID);

      // Sets the descriptor and unlocks the mutex.
      void FinalizeCreate(PHeapFileSystem::Descriptor* Desc);

      // Sets the descriptor and unlocks the mutex.
      // This version does not add reference counts for nested objects,
      // but adds RefCount-1 to this->ReferenceCount.  (RefCount-1 to allow for
      // the initial count from CreateLocked())
      void FinalizeCreate(PHeapFileSystem::Descriptor* Desc, int RefCount);

      template <typename T>
      T* Lock();

      void AddLock();
      bool SubLock();

      void AddReference();
      bool SubReference();

      //      PHeapObject* DeepCopy();

      // precondition: caller has a lock
      template <typename T>
      PHeapObject* CopyOnWrite(T*& Value);

      // precondition: object is locked.
      // removes the disk-copy of the object (if any).
      void SetDirty();

      id_type ID() const { return ObjectID; }

      // copies the object to the given stream - this function is only used once,
      // needs to be re-factored.  If we serialize the object, then write the format byte.
      // if we are copying disk descriptors, the format byte is already there.
      void Write(PHeapFileSystem::opheapstream& Out);

      // Adds Count to the reference count.  Count must be >= 0.
      void AddReference(int Count);

      // Flushes the object to disk (if it wasn't already), sets the internal
      // descriptor to NULL and returns the real descriptor.  
      // This has the effect of deleting the PHeapObject.
      PHeapFileSystem::Descriptor Persist();

      // for debug use, prints info on the internal contents to the specified stream.
      void DebugPrint(std::ostream& out) const;

      // forces removal of the descriptor (if it is non-NULL).  This results in
      // the object being in an invalid state, and further use is likely to be catastrophic.
      // Useful only when the underlying file system is about to be closed and the
      // object is not wanted.  If the lock count is non-zero, then it is incremented
      // (meaning the object in memory will never be destroyed).  Also prints a warning
      // message on std::cerr.
      void EmergencyDelete();

      // Internal use only, flushes the object to disk in the case that no file system was available
      // at the time that the lock count went to zero.
      void DoPendingFlush();

   private:
      PHeapObject(); // not implemented
      PHeapObject(PHeapObject const&); // not implemented
      PHeapObject& operator=(PHeapObject const&); // not implemented

      explicit PHeapObject(Private::PointerToObject* Obj, id_type ID);

      PHeapObject(PHeapFileSystem::Descriptor* Desc, id_type ID, 
		  int InitialReferenceCount);

      ~PHeapObject();

      bool DoLockCountZero(); // the slow part of SubLock()

      PHeapObject* CopyOnWriteLockCountZero();  // the non-template part of CopyOnWrite()

      AtomicRefCount ReferenceCount;
      AtomicRefCount LockCount;

      bool PendingFlush; // true if the lock count is zero but we cannot flush to disk
      // because there is no writable file system yet.

      id_type ObjectID; // this is set at creation and never modified

      // ObjectMutex protects access to Object and Heap.  Object is readable without
      // the mutex as long as the lock count is guaranteed to be > 0 for the duration of the access.
      pthread::mutex ObjectMutex;
      Private::PointerToObject* Object;   // non-NULL if the object is in memory
      PHeapFileSystem::Descriptor* MyDescriptor; // non-NULL if this object is stored on the heap, protected by ObjectMutex
};

inline
std::ostream& operator<<(std::ostream& out, PHeapObject const& Obj)
{
   Obj.DebugPrint(out);
   return out;
}

// Allocates a new object.  The initial state is lock count of 1.
template <typename T>
PHeapObject* Allocate(T* Ptr);

// if ID exists in the heap, then a reference is added and the object is returned.
// Otherwise returns NULL.
PHeapObject* AddReference(id_type ID);

// Shortcut for Obj = AddReference(ID); Ptr = Obj->Lock(); Obj->SubReference();
template <typename T>
std::pair<T*, PHeapObject*> GetObject(id_type ID)
{
   TRACE_PHEAP("GetObject()")(ID);
   if (ID == 0)
   {
      return std::pair<T*, PHeapObject*>(NULL, NULL);
   }

   PHeapObject* Obj = AddReference(ID);
   T* Ptr = Obj->template Lock<T>();
   Obj->SubReference();
   return std::pair<T*, PHeapObject*>(Ptr, Obj);
}

// Initializes the persistent heap
void Initialize(std::string const& FileName, int NumFiles, 
                size_t PageSize, size_t PageCacheByteSize,
                bool Unlink = false, bool AllowOverwrite = false);

// sets the format to use for persistent writing.  Defaults to PStream::format::Host
void SetPHeapFormat(int f);

// returns the format to be used for writing
int PHeapFormat();

// initializes the persistent heap from a given filesystem
PHeapObject* OpenPersistent(std::string const& FileName, size_t PageCacheByteSize, bool ReadOnly = false);

// only allow reading persisent files that have the specified page file metadata version.  Set to -1
// to enable all possible versions
void SetExpectedPageFileVersion(int v);

// returns the expected page file version.  -1 indicates no prefered version.
int ExpectedPageFileVersion();

// returns the page size of the current heap, or zero if there is no heap yet established.
size_t CurrentPageSize();

// Sets the page size to use as a default, for calls which require a page size parameter
void SetDefaultPageSize(size_t PSize);

// returns the default page size.  The initial value is 8 megabytes.
size_t DefaultPageSize();

// shuts down the heap without saving anything, deletes the file system.
void Shutdown();

// shuts down the persistent heap, flushing all data and writing the metadata.
void ShutdownPersistent(PHeapObject* MainObject);

// This is used to clean up files on termination.  If the persistent heap has
// already been deallocated then do nothing (so this is safe to call as the last function
// before main() exits).  If the persistent heap still exists and was created new, then
// delete the associated files (they won't be readable anyway if the heap hasn't been shutdown
// persistently).  If the persistent heap still exists but was loaded from an existing file,
// then do nothing - the file should still be readable after termination.
void Cleanup();

// exports an object to another filesystem, returning the PageId where the
// metadata for the heap is stored.  This function is useful only to
// implement higher level export functions.
PHeapFileSystem::PageId ExportHeap(PHeapFileSystem::BlockFileSystem* FS_, PHeapObject* MainObject);

// Exports an object to the new persistent heap of the given filename.
// The default PageSize is CurrentPageSize(), or
// DefaultPageSize() if CurrentPageSize() == 0.
void ExportHeap(std::string const& FileName, PHeapObject* MainObject, 
		int NumFiles = 1, size_t PageSize = 0);

// imports an object from another filesystem.
PHeapObject* ImportHeap(std::string const& FileName);

// returns the number of allocated pheap objects.  
// Useful diagnostic for detecting leaks.
size_t PHeapSize();

class Loader
{
   public:
      virtual ~Loader() {}

      // Loads object ID and returns a new'ed descriptor.
      // Nested objects should NOT have their reference count increased.
      virtual PHeapFileSystem::Descriptor* Load(id_type ID, PHeapFileSystem::BlockFileSystem* FS_) = 0;
};

PHeapObject* Inject(id_type ID, Loader* L);

void DebugHeap();

// helper function to split path files
inline
std::pair<std::string, std::string>
SplitPathFile(std::string const& PathFile)
{
   std::string::const_iterator PathSep = PathFile.begin();
   std::string::const_iterator Next = std::find(PathSep, PathFile.end(), '/');
   while (Next != PathFile.end())
   {
      ++Next;
      PathSep = Next;
      Next = std::find(PathSep, PathFile.end(), '/');
   }

   std::string Path = std::string(PathFile.begin(), PathSep);
   std::string File = std::string(PathSep, PathFile.end());
   return std::pair<std::string, std::string>(Path, File);
}

} // namespace pheap

#include "pheap.cc"

#endif
