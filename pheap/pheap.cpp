// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheap.cpp
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

#include "config.h"
#include "pheap.h"
#include "common/hash.h"
#include "common/atomicrefcount.h"
#include "pagestream.h"
#include <iomanip>
#include <set>
#include <atomic>

namespace pheap
{

using namespace PHeapFileSystem;

MessageLogger::Logger PHeapLog("PHEAPFS", std::cout);

namespace Private
{

typedef ext::hash_map<id_type, PHeapObject*> GlobalHeapType;

// we cannot control order of calling static constructors, so all of the important
// data for the persistent heap is instead stored on the heap and initialized
// by the nifty counter.  If objects are created prior to Initialize() being called,
// they are added to the PendingFlushList instead of being written to disk.  During
// Initialize(), the PendingFlushList is written to disk.
// The PendingFlushList is also used if an object wants to be flushed but
// the filesystem is read-only.

struct GlobalHeapDataType
{
   GlobalHeapType GlobalHeap;
   pthread::mutex HeapMutex;
   BlockFileSystem* FileSystem;
   std::atomic<id_type> SequenceNumber;
   id_type InitialSequenceNumber;
   pthread::mutex PendingFlushMutex;
   std::set<PHeapObject*> PendingFlushList;
   size_t MyDefaultPageSize;
   int DefaultFormat;
   std::list<PageId> MetadataPages;

   static GlobalHeapDataType* Data;
};

GlobalHeapDataType* GlobalHeapDataType::Data = NULL;

//
// GetInitialSequenceNumber
// Returns a pseudo-unique 64-bit number to use as the initial sequence number
// for object identifiers.
//
id_type GetInitialSequenceNumber()
{
   inttype::uint32 Unique = ext::get_unique_no_mutex();
   id_type isn = Unique;
   isn <<= 32;
   isn += ext::get_unique_no_mutex() & 0xFFF00000;
   // lower 20 bits are zero - does this reduce probability of collision?
   return isn;
}

void InitializeGlobalData()
{
   PRECONDITION(GlobalHeapDataType::Data == NULL);
   GlobalHeapDataType::Data = new GlobalHeapDataType();
   GlobalHeapDataType::Data->FileSystem = NULL;
   GlobalHeapDataType::Data->InitialSequenceNumber = GetInitialSequenceNumber();
   GlobalHeapDataType::Data->SequenceNumber.
     exchange(GlobalHeapDataType::Data->InitialSequenceNumber);
   GlobalHeapDataType::Data->MyDefaultPageSize = 8192 * 1024;  // 8 megabytes
   GlobalHeapDataType::Data->DefaultFormat = PStream::format::Host;
}

// returns true if the file system is initialized and not read-only.
bool IsFileSystemWritable()
{
   return !(GlobalHeapDataType::Data->FileSystem == NULL ||
            GlobalHeapDataType::Data->FileSystem->is_read_only());
}

BlockFileSystem*& FileSystem()
{
   return GlobalHeapDataType::Data->FileSystem;
}

GlobalHeapType& GlobalHeap()
{
   return GlobalHeapDataType::Data->GlobalHeap;
}

pthread::mutex& GlobalHeapMutex()
{
   return GlobalHeapDataType::Data->HeapMutex;
}

std::atomic<id_type>& SequenceNumber()
{
  return GlobalHeapDataType::Data->SequenceNumber;
}

size_t& GlobalDefaultPageSize()
{
   return GlobalHeapDataType::Data->MyDefaultPageSize;
}

std::list<PageId>&
GlobalMetadataPages()
{
   return GlobalHeapDataType::Data->MetadataPages;
}

PHeapObject* PendingFlushListTop()
{
   pthread::mutex::sentry(GlobalHeapDataType::Data->PendingFlushMutex);
   return GlobalHeapDataType::Data->PendingFlushList.empty() ? NULL :
    *GlobalHeapDataType::Data->PendingFlushList.begin() ;
}

void SetPendingFlush(PHeapObject* Obj)
{
   pthread::mutex::sentry(GlobalHeapDataType::Data->PendingFlushMutex);
   GlobalHeapDataType::Data->PendingFlushList.insert(Obj);
}

void ClearPendingFlush(PHeapObject* Obj)
{
   pthread::mutex::sentry(GlobalHeapDataType::Data->PendingFlushMutex);
   GlobalHeapDataType::Data->PendingFlushList.erase(Obj);
}

// utility function to make a deep copy of a Descriptor into another buffer.
void DeepCopyDescriptor(PStream::generic_opstreambuf* Buf, Descriptor const& Desc)
{
   // copy the data part across
   Descriptor::block_iterator CurrentBlock = Desc.block_begin();
   while (CurrentBlock != Desc.block_end())
   {
      InputBuffer InBuf(*CurrentBlock++);
      Buf->put_n(InBuf.buf_ptr(), InBuf.remain());
   }

   // copy the object id's across
   Descriptor::id_iterator CurrentId = Desc.id_begin();
   while (CurrentId != Desc.id_end())
   {
      Buf->put_id(*CurrentId++);
   }
}

// utility function for copying a Descriptor to a different file system
Descriptor CopyDescriptorToFileSystem(Descriptor const& Desc, BlockFileSystem* FS_)
{
   Descriptor Result;
   opheapstream Buf(FS_, &Result);  // ** No format here, that is already in the descriptor
   DeepCopyDescriptor(Buf.get_buffer(), Desc);
   return Result;
}

// HeapType is a hash_map of id_type and HeapRecord.  This
// represents the heap metadata.  It is directly streamable.
struct HeapRecord
{
   HeapRecord() : RefCount(0) {}
   int RefCount;
   Descriptor Desc;
};

void WriteHeapRecord(BlockFileSystem* FS_, PStream::opstream& out, HeapRecord const& Rec)
{
   out.write<uint32>(Rec.RefCount);
   Rec.Desc.Write(FS_, out);
}

void ReadHeapRecord(BlockFileSystem* FS_, PStream::ipstream& in, HeapRecord& Rec)
{
   Rec.RefCount = in.read<uint32>();
   Rec.Desc = Descriptor(FS_, in);
}

typedef ext::hash_map<id_type, HeapRecord> HeapType;  // the metadata container

void WriteHeap(BlockFileSystem* FS_, PStream::opstream& out, HeapType const& Heap)
{
   PRECONDITION(FS_ != NULL);
   out.write<uint32>(Heap.size());
   for (HeapType::const_iterator I = Heap.begin(); I != Heap.end(); ++I)
   {
      out.write<uint64>(I->first);
      WriteHeapRecord(FS_, out, I->second);
   }
}

void ReadHeap(BlockFileSystem* FS_, PStream::ipstream& in, HeapType& Heap)
{
   PRECONDITION(FS_ != NULL);
   Heap.clear();
   size_t Size = in.read<uint32>();
   for (size_t i = 0; i < Size; ++i)
   {
      id_type ID = in.read<uint64>();
      HeapRecord Rec;
      ReadHeapRecord(FS_, in, Rec);
      Heap[ID] = Rec;
   }
}

id_type AllocateID()
{
   return SequenceNumber()++;
}

BlockFileSystem* GetFileSystem() { return GlobalHeapDataType::Data->FileSystem; }

// This is the initialization function called by the nifty counter
void InitializePHeap()
{
   InitializeGlobalData();
}

} // namespace Private

using namespace Private;

//
// PHeapObject
//

int PHeapFormat()
{
   return Private::GlobalHeapDataType::Data->DefaultFormat;
}

void SetPHeapFormat(int f)
{
   Private::GlobalHeapDataType::Data->DefaultFormat = f;
}

// add the references for nested objects.
void AddNestedReferences(Descriptor* Desc)
{
   PRECONDITION(Desc != NULL);
   // we can do everything with the GlobalHeapMutex
   pthread::mutex::sentry Lock(GlobalHeapMutex());
   for (Descriptor::id_iterator I = Desc->id_begin(); I != Desc->id_end(); ++I)
   {
      if (*I != 0)
      {
         PHeapObject* Obj = GlobalHeap()[*I];
         CHECK(Obj != NULL);
         Obj->AddReference();
      }
   }
}

// remove the references for nested objects.
void SubNestedReferences(Descriptor* Desc)
{
   PRECONDITION(Desc != NULL);
   if (!Desc) return;
   PHeapObject* Obj;
   // This is a bit cumbersome as we have to grab the GlobalHeapMutex to lookup the
   // GlobalHeap, but we have to give up the mutex before calling SubReference()
   for (Descriptor::id_iterator I = Desc->id_begin(); I != Desc->id_end(); ++I)
   {
      if (*I != 0)
      {
         {
            pthread::mutex::sentry Lock(GlobalHeapMutex());
            Obj = GlobalHeap()[*I];
         }
         CHECK(Obj != NULL);
         Obj->SubReference();
      }
   }
}

PHeapObject::~PHeapObject()
{
   TRACE_PHEAP("Deleting ")(this)(MyDescriptor);
   delete Object;
   if (MyDescriptor)
   {
      SubNestedReferences(MyDescriptor);
      TRACE_PHEAP(MyDescriptor)(*MyDescriptor);
      delete MyDescriptor;
   }
   pthread::mutex::sentry Lock(GlobalHeapMutex());
   GlobalHeap().erase(ObjectID);
}

PHeapObject::PHeapObject(Private::PointerToObject* Obj, id_type ID)
   : ReferenceCount(1), LockCount(1), PendingFlush(false), ObjectID(ID), Object(Obj), MyDescriptor(NULL)
{
   CHECK(ID != 0);
}

PHeapObject::PHeapObject(Descriptor* Desc, id_type ID, int InitialReferenceCount)
  : ReferenceCount(InitialReferenceCount), LockCount(0), PendingFlush(false),
    ObjectID(ID), Object(NULL), MyDescriptor(Desc)
{
   CHECK(ID != 0);
}

PHeapObject* PHeapObject::Create(Private::PointerToObject* Obj, id_type ID)
{
   pthread::mutex::sentry Lock(GlobalHeapMutex());
   CHECK(GlobalHeap().count(ID) == 0);
   return GlobalHeap()[ID] = new PHeapObject(Obj, ID);
}

PHeapObject* PHeapObject::Create(Descriptor const& Desc, id_type ID, int InitialReferenceCount)
{
   pthread::mutex::sentry Lock(GlobalHeapMutex());
   PHeapObject*& Obj = GlobalHeap()[ID];
   CHECK(Obj == NULL);
   Obj = new PHeapObject(new Descriptor(Desc), ID, InitialReferenceCount);
   return Obj;
}

PHeapObject* PHeapObject::CreateLocked(id_type ID)
{
   CHECK(GlobalHeap()[ID] == NULL);
   PHeapObject* Obj = new PHeapObject(NULL, ID, 1);
   GlobalHeap()[ID] = Obj;
   Obj->ObjectMutex.lock();
   return Obj;
}

void PHeapObject::FinalizeCreate(Descriptor* Desc)
{
   CHECK(MyDescriptor == NULL);
   CHECK(Object == NULL);
   MyDescriptor = Desc;
   AddNestedReferences(MyDescriptor);
   ObjectMutex.unlock();
}

void PHeapObject::FinalizeCreate(Descriptor* Desc, int RefCount)
{
   CHECK(MyDescriptor == NULL);
   CHECK(Object == NULL);
   MyDescriptor = Desc;
   AddReference(RefCount-1);
   ObjectMutex.unlock();
}

bool PHeapObject::DoLockCountZero()
{
   TRACE_PHEAP("PHeapObject::DoLockCountZero()")(*this);
   //ObjectMutex.lock(); // ** moved to SubLock function

   if (LockCount.is_zero())
   {
      // we are going to destroy the memory object for sure,
      // so we decrement the reference count too.
      if (--ReferenceCount == 0)
      {
         //ObjectMutex.unlock();
         delete this;
         return true;
      }
      else
      {
         // If we don't have a descriptor then the object is not currently on disk
         if (!MyDescriptor)
         {
            TRACE_PHEAP("PHeapObject::DoLockCountZero(): no descriptor - writing to disk");
            //      DEBUG_TRACE("MyDescriptor is NULL");
            // if we can't write to the file system, then add this to the set of pending writes.
            if (!IsFileSystemWritable())
            {
               TRACE_PHEAP("PHeapObject::DoLockCountZero(): Filesystem is not writable")(*this);
               PendingFlush = true;
               SetPendingFlush(this);
               //ObjectMutex.unlock();
               return false; // NOTE: we don't delete the Object here
            }
            MyDescriptor = new Descriptor();
            {
               opheapstream Out(FileSystem(), MyDescriptor); Out.put_format(PHeapFormat());
               Object->WriteToStream(Out);
            }
            AddNestedReferences(MyDescriptor);
         }
         delete Object; Object = NULL;
      }
   }
   //ObjectMutex.unlock();
   return false;
}

void PHeapObject::AddReference(int Count)
{
   for (int i = 0; i < Count; ++i)
   {
      ++ReferenceCount;
   }
}

PHeapObject* PHeapObject::CopyOnWriteLockCountZero()
{
   // the lock count is 1, we are the only user of Object.
   // Are we the only holder of a reference too?
   if (--ReferenceCount == 0)
   {
      // yes, we can re-use this PHeapObject.  We still need to change the ID though.
      pthread::mutex::sentry Lock(GlobalHeapMutex());
      GlobalHeap().erase(ObjectID);
      //      std::cout << " old ID " << ObjectID;
      ObjectID = Private::AllocateID();
      //      std::cout << " new ID " << ObjectID << std::endl;
      GlobalHeap()[ObjectID] = this;
      ++ReferenceCount;
      ++LockCount;
      return this;
   }
   else
   {
      // we need to move Object to a new PHeapObject, and leave a disk-copy.
      if (!MyDescriptor)
      {
         // if the FileSystem is not writable, then we cannot reuse the memory object, but instead
         // we make a deep copy of it and leave the old version intact, and add it to the list of
         // pending objects.
         if (!IsFileSystemWritable())
         {
            PendingFlush = true;
            SetPendingFlush(this);
            PHeapObject* Other = Create(Object->clone(), Private::AllocateID());

            // The swap here is needed because whoever called CopyOnWrite might have a copy of
            // our Object pointer somewhere, so we shouldn't change it under them.  So,
            // return a new PHeapObject that contains the old Object pointer.
            // The swap should also be safe, because the only other things pointing to this
            // PHeapObject are references, not locks, and therefore should not have a copy of the Object pointer.
            std::swap(Object, Other->Object);

            return Other;
         }
         MyDescriptor = new Descriptor();
         {
            opheapstream Out(FileSystem(), MyDescriptor); Out.put_format(PHeapFormat());
            Object->WriteToStream(Out);
         }
         AddNestedReferences(MyDescriptor);
      }
      PHeapObject* Other = Create(Object, Private::AllocateID());
      Object = NULL;
      return Other;
   }
}

void PHeapObject::Write(PHeapFileSystem::opheapstream& Out)
{
   pthread::mutex::sentry Lock(ObjectMutex);

   if (Object)
   {
      Out.put_format(PHeapFormat());
      Object->WriteToStream(Out);
   }
   else if (MyDescriptor)
   {
      DeepCopyDescriptor(Out.get_buffer(), *MyDescriptor);
   }
   else
   {
      PANIC("Peristent object was neither in memory or on the heap!");
   }
}

void PHeapObject::DoPendingFlush()
{
   ObjectMutex.lock();
   CHECK(PendingFlush);
   CHECK(MyDescriptor == NULL);
   CHECK(Object != NULL);
   MyDescriptor = new Descriptor();
   {
      opheapstream Out(FileSystem(), MyDescriptor); Out.put_format(PHeapFormat());
      Object->WriteToStream(Out);
   }
   AddNestedReferences(MyDescriptor);
   delete Object; Object = NULL;
   ClearPendingFlush(this);
   PendingFlush = false;
   ObjectMutex.unlock();
   if (ReferenceCount.is_zero()) delete this;
}

Descriptor PHeapObject::Persist()
{
   Descriptor Desc;
   {
      pthread::mutex::sentry Lock(ObjectMutex);
      if (!MyDescriptor)
      {
         CHECK(Object != NULL);
         //      std::cerr << "WARNING: Object " << ObjectID << " type " << Object->GetTypeid().name()
         //                << " has memory references that are now invalid.\n";
         ++LockCount;  // prohibit deleting the object, trying to write it to disk would be very bad
         Descriptor Desc;
         {
            opheapstream Out(FileSystem(), &Desc); Out.put_format(PHeapFormat());
            Object->WriteToStream(Out);
         }
         return Desc;
      }
      ++ReferenceCount;   // stop the object being deallocated
      if (!LockCount.is_zero()) ++LockCount;
      Desc = *MyDescriptor;
      delete MyDescriptor;
      MyDescriptor = NULL;  // prevent the destructor from removing the file references
   }
   //   delete this;
   return Desc;
}

void PHeapObject::DebugPrint(std::ostream& out) const
{
   out << "  this: " << this << "  locked: " << LockCount.value()
       << "  ref: " << ReferenceCount.value()
       << "  object: " << (void*) Object;
   if (Object) out << "  type: " << Object->GetTypeid().name();
   out << "  Descriptor: " << (void*) MyDescriptor;
   if (MyDescriptor)
   {
      out << "  size: " << MyDescriptor->data_size() << "  nested: " << MyDescriptor->id_size();
      out << " detail:" << *MyDescriptor;
   }
}

void PHeapObject::SetDirty()
{
   pthread::mutex::sentry Lock(ObjectMutex);
   DEBUG_PRECONDITION(Object != NULL);
   if (MyDescriptor)
   {
      SubNestedReferences(MyDescriptor);
      delete MyDescriptor;
      MyDescriptor = NULL;
   }
}

void PHeapObject::EmergencyDelete()
{
   pthread::mutex::sentry Lock(ObjectMutex);
   //   std::cerr << "WARNING: Object " << ObjectID
   //        << " has memory references but no persistent references; the memory references are now invalid.\n";
   if (MyDescriptor)
   {
      delete MyDescriptor; MyDescriptor = NULL;
   }
   if (!LockCount.is_zero()) ++LockCount;
   ++ReferenceCount;
}

//
// globals
//

PHeapObject* AddReference(id_type ID)
{
   {
      pthread::mutex::sentry Lock(GlobalHeapMutex());
      if (GlobalHeap().count(ID) != 0)
      {
         PHeapObject* Obj = GlobalHeap()[ID];
         Obj->AddReference();
         return Obj;
      }
   }
   return NULL;
}

void Initialize(std::string const& FileName, int NumFiles, size_t PageSize,
                size_t PageCacheByteSize, bool Unlink, bool AllowOverwrite)
{
   pthread::mutex::sentry Lock(GlobalHeapMutex());

   PRECONDITION(FileSystem() == NULL);
   notify_log(30, PHeapLog) << "Initializing persistent storage.\n";
   FileSystem() = new PHeapFileSystem::BlockFileSystem();
   FileSystem()->create(FileName, NumFiles, PageSize, PageCacheByteSize, Unlink, AllowOverwrite);
   notify_log(40, PHeapLog) << "Initial sequence number is "
                         << GlobalHeapDataType::Data->InitialSequenceNumber << '\n';
   // walk the pending flush list.
   // This would mean that the object is already created and the lock count is zero,
   // so we can write it to disk now.
   PHeapObject* Obj;
   while ((Obj = PendingFlushListTop()) != NULL)
   {
      Obj->DoPendingFlush();
   }
}

PHeapObject* OpenPersistent(std::string const& FileName, size_t PageCacheByteSize, bool IsReadOnly)
{
   PRECONDITION(GetFileSystem() == NULL);
   FileSystem() = new PHeapFileSystem::BlockFileSystem();
   PageId MetaId = FileSystem()->open(FileName, PageCacheByteSize, IsReadOnly);

   ipagestream MetaIn(FileSystem(), MetaId, PStream::format::XDR);
   id_type MainObjectID = MetaIn.read<uint64>();

   HeapType HeapRecords;
   ReadHeap(FileSystem(), MetaIn, HeapRecords);

   // Save the metadata pages so that we don't overwrite them
   GlobalMetadataPages() = MetaIn.pages();

   for (HeapType::const_iterator I = HeapRecords.begin(); I != HeapRecords.end(); ++I)
   {
      PHeapObject::Create(I->second.Desc, I->first, I->second.RefCount);
   }

   pthread::mutex::sentry Lock(GlobalHeapMutex());
   return GlobalHeap()[MainObjectID];  // no need to increment the ref count, this is already done
}

void ShutdownPersistent(PHeapObject* MainObject)
{
   PRECONDITION(MainObject != NULL);
   //   pthread::mutex::sentry Lock(GlobalHeapMutex());
   PRECONDITION(GetFileSystem() != NULL);

   // we can delete the old metadata
   GlobalMetadataPages().clear();

   //   FileSystem()->defragment();

   // write the metadata

   HeapType HeapRecords;
   std::list<id_type> MissingIds;

   id_type MainObjectID = MainObject->ID();

   // set the main object as a missing object, with a reference count of 1
   HeapRecords[MainObject->ID()].RefCount = 1;
   MissingIds.push_back(MainObject->ID());

  typedef Descriptor::id_iterator IDIter;

   while (!MissingIds.empty())
   {
      id_type ID = MissingIds.front();
      MissingIds.pop_front();

      Descriptor Desc = HeapRecords[ID].Desc = GlobalHeap()[ID]->Persist();

      // increment the reference count for the nested objects, and
      // if any were previously not present (ref count of zero), add
      // the ID to the MissingIds list.
      for (IDIter I = Desc.id_begin(); I != Desc.id_end(); ++I)
      {
         if (*I != 0 && HeapRecords[*I].RefCount++ == 0) MissingIds.push_back(*I);
      }
   }

   // Look for any objects in the global heap that are not going to be saved
   for (GlobalHeapType::const_iterator I = GlobalHeap().begin(); I != GlobalHeap().end(); ++I)
   {
      if (HeapRecords.count(I->first) == 0)
      {
         I->second->EmergencyDelete();
      }
   }

   opagestream MetaOut(FileSystem(), PStream::format::XDR);
   MetaOut.write<uint64>(MainObjectID);
   WriteHeap(FileSystem(), MetaOut, HeapRecords);

   PageId MetaId = MetaOut.commit();
   FileSystem()->persistent_shutdown(MetaId);

   delete FileSystem();
   FileSystem() = NULL;
}

void Shutdown()
{
   pthread::mutex::sentry Lock(GlobalHeapMutex());
   PRECONDITION(GetFileSystem() != NULL);
   for (GlobalHeapType::const_iterator I = GlobalHeap().begin(); I != GlobalHeap().end(); ++I)
   {
      I->second->EmergencyDelete();
   }
   GlobalHeap().clear();
   FileSystem()->shutdown(true);  // remove the file
   delete FileSystem();
   FileSystem() = NULL;
}

void Cleanup()
{
   pthread::mutex::sentry Lock(GlobalHeapMutex());
   if (GetFileSystem() == NULL)
      return;

   FileSystem()->cleanup();
}

PageId ExportHeap(BlockFileSystem* FS_, PHeapObject* MainObject)
{
   PRECONDITION(FS_ != NULL);

   HeapType HeapRecords;

   // Keep a list of Id's that we need to add to the exported heap
   std::list<id_type> MissingIds;

   typedef Descriptor::id_iterator IDIter;

   // seed the list of id's to add, with the main object.  The reference count starts at 1,
   // to account for the main object pointer.
   HeapRecords[MainObject->ID()].RefCount = 1;
   MissingIds.push_back(MainObject->ID());

   // Inject the front of the MissingIds list into the heap, until
   // MissingIds is empty.
   while (!MissingIds.empty())
   {
      id_type ID = MissingIds.front();
      MissingIds.pop_front();

      // Get the object at ID and write it to the exported heap.
      PHeapObject* Obj = AddReference(ID);
      Descriptor* Desc = &HeapRecords[ID].Desc;
      opheapstream Out(FS_, Desc); //Out.put_format(PHeapFormat());
      Obj->Write(Out);
      Obj->SubReference();  // since we had to add a reference to get the PHeapObject*

      // increment the reference count for the nested objects, and
      // if any were previously not present (ref count of zero), add
      // the ID to the MissingIds list.
      for (IDIter I = Desc->id_begin(); I != Desc->id_end(); ++I)
      {
         if (*I != 0 && HeapRecords[*I].RefCount++ == 0) MissingIds.push_back(*I);
      }
   }

   opagestream MetaOut(FS_, PStream::format::XDR);
   MetaOut.write<uint64>(MainObject->ID());
   WriteHeap(FS_, MetaOut, HeapRecords);
   PageId MetaId = MetaOut.commit();
   HeapRecords.clear();
   return MetaId;
}

void ExportHeap(std::string const& FileName, PHeapObject* MainObject,
                int NumFiles, size_t PageSize)
{
   BlockFileSystem* FS_ = new BlockFileSystem;
   if (PageSize == 0)
   {
      PageSize = pheap::CurrentPageSize();
      if (PageSize == 0) PageSize = pheap::DefaultPageSize();
   }
   FS_->create(FileName, NumFiles, PageSize, PageSize*2);
   PageId MetaId = ExportHeap(FS_, MainObject);
   FS_->persistent_shutdown(MetaId);
   delete FS_;
}

PHeapObject* ImportHeap(BlockFileSystem* FS_, PageId MetaPage)
{
   PRECONDITION(FileSystem() != NULL);
   id_type MainObjectID;
   HeapType HeapRecords;

   ipagestream MetaIn(FS_, MetaPage, PStream::format::XDR);

   MainObjectID = MetaIn.read<uint64>();
   ReadHeap(FS_, MetaIn, HeapRecords);

   CHECK(MainObjectID != 0);

   GlobalHeapMutex().lock();

   for (HeapType::const_iterator I = HeapRecords.begin(); I != HeapRecords.end(); ++I)
   {
      // if the object is not present, allocate it.
      if (GlobalHeap().count(I->first) == 0)
      {
         // copy the Descriptor into MyFileSystem.  Do it as a two-stage
         // CreateLocked() / FinalizeCreate() so that we can drop the GlobalHeapMutex
         // while copying the Descriptor to the new filesystem.  This should
         // allow multiple concurrent calls to ImportHeap() to be efficient.
         PHeapObject* Obj = PHeapObject::CreateLocked(I->first);
         GlobalHeapMutex().unlock();
         Descriptor* Desc = new Descriptor(CopyDescriptorToFileSystem(I->second.Desc, FileSystem()));
         Obj->FinalizeCreate(Desc, I->second.RefCount);
         GlobalHeapMutex().lock();
      }
      else
      {
         // otherwise, just bump the reference count.
         // It might be desirable to overwrite the previous version, but
         // We can't really do anything else as the object may be locked.
         CHECK(GlobalHeap()[I->first] != NULL);
         GlobalHeap()[I->first]->AddReference(I->second.RefCount);
      }
   }
   GlobalHeapMutex().unlock();

   HeapRecords.clear();  // this deletes the Descriptors so we can delete FS_ properly.

   pthread::mutex::sentry Lock(GlobalHeapMutex());
   return GlobalHeap()[MainObjectID]; // no need to add a reference, this is done implicitly already
}

// This is the version that is defined in the header
PHeapObject* ImportHeap(std::string const& File)
{
   BlockFileSystem* NewFS = new BlockFileSystem;
   PHeapObject* Obj = ImportHeap(NewFS, NewFS->open(File, 0, true));
   NewFS->shutdown();
   delete NewFS;
   return Obj;
}

PHeapObject* Inject(id_type ID, Loader* L)
{
   PRECONDITION(FileSystem() != NULL);
   // There are a few ways to do this.  A safe way is to grab the GlobalHeapMutex and hold onto
   // it until all nested objects have been constructed.  This is a bit rude in a multithread
   // environment though, as another thread that called, say, AddReference(ID) for an unrelated
   // object would block.  Instead, we create the PHeapObject such that it initially has its
   // mutex locked.  This means that we don't need to hold the GlobalHeapMutex, and the only
   // operations what will block for a length of time are operations that require the loaded
   // object itself.

   // See if the object is already present.  If it is, just add a reference.
   // If it isn't, we need to create it and load it.

   PHeapObject* Obj;
   {
      pthread::mutex::sentry Lock(GlobalHeapMutex());
      Obj = GlobalHeap()[ID];
      if (Obj)
      {
         Obj->AddReference();
         return Obj;
      }
      // else
      Obj = PHeapObject::CreateLocked(ID);
   }

   Descriptor* Desc = L->Load(ID, FileSystem());
   for (Descriptor::id_iterator I = Desc->id_begin(); I != Desc->id_end(); ++I)
   {
      if (*I != 0) Inject(*I, L);
   }
   Obj->FinalizeCreate(Desc, 1); // don't add nested references, we've already done it.
   return Obj;
}

void DebugHeap()
{
   std::cerr << "\npheap debug report:\n"
       << "Number of objects: " << GlobalHeap().size() << "\n";
   for (GlobalHeapType::const_iterator I = GlobalHeap().begin(); I != GlobalHeap().end(); ++I)
   {
     std::cerr << " identifier: " << std::setw(20) << I->first << "  hashkey: "
         << std::setw(5) << GlobalHeap().hash_funct()(I->first)
         << "  pointer: " << (void*) I->second;
      I->second->DebugPrint(std::cerr);
      std::cerr << '\n';
   }
   FileSystem()->Debug();
   std::cerr << "end pheap debug report." << std::endl;
}

size_t PHeapSize()
{
   pthread::mutex::sentry Lock(GlobalHeapMutex());
   return GlobalHeap().size();
}

size_t CurrentPageSize()
{
   return GetFileSystem() ? GetFileSystem()->get_page_size() : 0;
}

void SetDefaultPageSize(size_t PSize)
{
   // round PSize up to the next multiple of BufferAllocator::PageGranularity
   size_t Margin = PSize % BufferAllocator::PageGranularity;
   if (Margin != 0) PSize = PSize + BufferAllocator::PageGranularity - Margin;
   GlobalDefaultPageSize() = PSize;
}

size_t DefaultPageSize()
{
   return GlobalDefaultPageSize();
}

} // namespace pheap
