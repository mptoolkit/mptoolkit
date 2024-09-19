// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pheap/pheap.cc
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
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

#include "pheapstream.h"
#include "pheapallocator.h"

namespace pheap
{

void SubNestedReferences(PHeapFileSystem::Descriptor* Desc);

namespace Private
{

id_type AllocateID();

BlockFileSystem* GetFileSystem();

void ClearPendingFlush(PHeapObject* Obj);

//
// concrete classes derived from PointerToObject
//

template <typename T>
class PointerToValue : public PointerToObject
{
   public:
     explicit PointerToValue(T* Ptr_);

      virtual ~PointerToValue();
      virtual PointerToObject* clone() const;

      virtual void WriteToStream(PStream::opstream& out) const;

      virtual std::type_info const& GetTypeid() const;

      T* Get() const { return Ptr; }

      //      void* operator new(size_t size);
      //      void operator delete(void* p);

   private:
      T* Ptr;
      PointerToValue();  // not defined
};

template <typename T>
PointerToValue<T>::PointerToValue(T* Ptr_)
  : Ptr(Ptr_)
{
}

template <typename T>
PointerToValue<T>::~PointerToValue()
{
   delete Ptr;
}

template <typename T>
PointerToObject* PointerToValue<T>::clone() const
{
   return new PointerToValue<T>(new T(*Ptr));
}

template <typename T>
void PointerToValue<T>::WriteToStream(PStream::opstream& out) const
{
   out << *Ptr;
}

template <typename T>
std::type_info const& PointerToValue<T>::GetTypeid() const
{
   return typeid(T);
}

#if 0
template <typename T>
void* PointerToValue<T>::operator new(size_t size);
{
   DEBUG_PRECONDITION(size == sizeof(PointerToValue<T>));
   return PoolAlloc::allocate(size);
}

template <typename T>
void  PointerToValue<T,>::operator delete(void* p)
{
   PoolAlloc::deallocate(p, sizeof(PointerToValue<T, false>));
}
#endif

template <class BaseClass>
class PointerToPolymorphic : public PointerToObject
{
   public:
      explicit PointerToPolymorphic(BaseClass* Ptr_);

      virtual ~PointerToPolymorphic();
      virtual PointerToObject* clone() const;

      virtual void WriteToStream(PStream::opstream& out) const;

      virtual std::type_info const& GetTypeid() const;

      BaseClass* Get() const { return Ptr; }

      //      void* operator new(size_t size);
      //      void operator delete(void* p);

   private:
      BaseClass* Ptr;
      PointerToPolymorphic();  // not defined
};

template <class BaseClass>
PointerToPolymorphic<BaseClass>::PointerToPolymorphic(BaseClass* Ptr_)
  : Ptr(Ptr_)
{
}

template <class BaseClass>
PointerToPolymorphic<BaseClass>::~PointerToPolymorphic()
{
   delete Ptr;
}

template <class BaseClass>
PointerToObject* PointerToPolymorphic<BaseClass>::clone() const
{
   return new PointerToPolymorphic<BaseClass>(Ptr->clone());
}

template <class BaseClass>
 void PointerToPolymorphic<BaseClass>::WriteToStream(PStream::opstream& out) const
{
   out << *Ptr;
}

template <class BaseClass>
std::type_info const& PointerToPolymorphic<BaseClass>::GetTypeid() const
{
   return typeid(*Ptr);
}

#if 0
template <class BaseClass>
void* PointerToPolymorphic<BaseClass>::operator new(size_t size)
{
   DEBUG_PRECONDITION(size == sizeof(PointerToPolymorphic<BaseClass>));
   return PoolAlloc::allocate(size);
}

template <class BaseClass>
void PointerToPolymorphic<BaselCass>::operator delete(void* p)
{
   PoolAlloc::deallocate(p, sizeof(PointerToPolymorphic<BaseClass>));
}
#endif

// a helper for PointerToObject

template <typename T, bool IsPolymorphic = PStream::poly_traits<T>::is_polymorphic>
struct ObjectHelper;


template <typename T>
struct ObjectHelper<T, true>
{
   static T* TryCastPointer(PointerToObject const* Obj);
   static PointerToObject* Create(T* Obj);
};

template <typename T>
struct ObjectHelper<T, false>
{
   static T* TryCastPointer(PointerToObject const* Obj);
   static PointerToObject* Create(T* Obj);
};

template <typename T>
T* ObjectHelper<T, false>::TryCastPointer(PointerToObject const* Obj)
{
   PointerToValue<T> const* AsValue = dynamic_cast<PointerToValue<T> const*>(Obj);
   CHECK(AsValue != NULL);
   CHECK(AsValue->Get() != NULL);
   return AsValue->Get();
}

template <typename T>
inline
PointerToObject* ObjectHelper<T, false>::Create(T* Obj)
{
   return new PointerToValue<T>(Obj);
}

template <class T>
T* ObjectHelper<T, true>::TryCastPointer(PointerToObject const* Obj)
{
   typedef typename T::heirachy_base BaseClass;
   PointerToPolymorphic<BaseClass> const* AsPoly =
      dynamic_cast<PointerToPolymorphic<BaseClass> const*>(Obj);
   T* Ptr = dynamic_cast<T*>(AsPoly->Get());
   CHECK(Ptr != NULL);
   return Ptr;
}

template <typename T>
inline
PointerToObject* ObjectHelper<T, true>::Create(T* Obj)
{
   typedef typename T::heirachy_base BaseClass;
   return new PointerToPolymorphic<BaseClass>(Obj);
}

//
// PointerToObject
//

template <typename T>
inline
T* PointerToObject::GetPointer() const
{
   return Private::ObjectHelper<T>::TryCastPointer(this);
}

template <typename T>
inline
PointerToObject* PointerToObject::Create(T* Obj)
{
   return Private::ObjectHelper<T>::Create(Obj);
}

} // namespace Private

//
// PHeapObject
//

template <typename T>
inline
PHeapObject* PHeapObject::Create(T* Ptr)
{
   return Create(Private::PointerToObject::Create(Ptr), Private::AllocateID());
}

template <typename T>
inline
PHeapObject* PHeapObject::Create(T* Ptr, id_type ID)
{
   return Create(Private::PointerToObject::Create(Ptr), ID);
}

template <typename T>
T* PHeapObject::Lock()
{
   // possible race condition: two threads call Lock() at the same time.
   // Averted by grabbing the mutex before anything happens.

   // possible race condition: One thread calls Lock() immediately after
   // another thread destroys a PHeapPointer which reduces the lock count to zero
   // but the object is still in memory.
   // Averted by the other thread re-checking the LockCount after grabbing the mutex.

   TRACE_PHEAP("PHeapObject::Lock()")(*this);

   std::lock_guard<std::mutex> Lock(ObjectMutex);

   // we cannot simply increment the lock count here,
   // we need to make sure that the object is in memory first.
   if (LockCount.is_zero())
   {
      // an object in memory counts as a reference.  We increment regardless of PendingFlush status
      ++ReferenceCount;
      if (Object == NULL)
      {
         // load the object into memory.
         TRACE_PHEAP("Loading object")(this);
         PHeapFileSystem::ipheapstream In(Private::GetFileSystem(), MyDescriptor);
         T* Ptr = PStream::ConstructFromStream<T>(In);
         Object = Private::PointerToObject::Create(Ptr);

         // If the filesystem is not writable, then we might as well delete the disk image
         if (Private::GetFileSystem()->is_read_only())
         {
            //SubNestedReferences(MyDescriptor);
            //TRACE_PHEAP(MyDescriptor);
            //delete MyDescriptor;
         }

         // we cannot increment the lock count until the object is loaded.
         // we also need a memory barrier here so that the object's memory
         // is coherent before the lock count becomes non-zero.
         //         write_memory_barrier();
         ++LockCount;

         return Ptr;
      }

      // if we get here, then the lock count is zero but the object is in memory.  The most likely
      // reason is a pending flush.
      if (PendingFlush && !Private::GetFileSystem()->is_read_only())
      {
         Private::ClearPendingFlush(this);
         PendingFlush = false;
      }
   }

   // if we get here then the object is already in memory.
   ++LockCount;
   return Object->template GetPointer<T>();
}

inline
void PHeapObject::AddReference()
{
   TRACE_PHEAP_X("PHeapObject::AddReference()")(*this);
   ++ReferenceCount;
}

inline
bool PHeapObject::SubReference()
{
   TRACE_PHEAP_X("PHeapObject::SubRference()")(*this);
   if (--ReferenceCount == 0)
   {
      delete this;
      return true;
   }
   return false;
}

inline
void PHeapObject::AddLock()
{
   TRACE_PHEAP_X("PHeapObject::AddLock()")(*this);
   ++LockCount;
}

inline
bool PHeapObject::SubLock()
{
   TRACE_PHEAP_X("PHeapObject::SubLock()")(*this);
   bool Result = false;
   ObjectMutex.lock(); // we have to grab the lock here because it isn't safe for two threads to enter DoLockCountZero()
   if (--LockCount == 0)
   {
      Result = this->DoLockCountZero();
   }
   ObjectMutex.unlock();
   if (Result)
      delete this;
   return Result;
}

template <typename T>
PHeapObject* PHeapObject::CopyOnWrite(T*& Value)
{
   std::lock_guard<std::mutex> Lock(ObjectMutex);
   DEBUG_CHECK(Object != NULL && Object->template GetPointer<T>() == Value);

   if (--LockCount == 0)
   {
      return CopyOnWriteLockCountZero();
   }
   else
   {
      // we are not the only one with a lock.  We have to copy the object itself.
      PHeapObject* Other = Create(Object->clone(), Private::AllocateID());
      Value = Other->Object->template GetPointer<T>();
      return Other;
   }
}

//
// namespace-level free functions
//

template <typename T>
inline
PHeapObject* Allocate(T* Ptr)
{
   return PHeapObject::Create(Ptr, Private::AllocateID());
}

template <typename T>
std::pair<T*, PHeapObject*> GetObject(id_type ID);

} // namespace pheap
