// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/prefptr.h
//
// Copyright (C) 1999-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  Created 1999-07-14
  Ian McCulloch

  Defines the pref_ptr and const_pref_ptr template classes.  These are reference counted
  pointers that exist on the normal heap, with capability to be written out to the persistent heap.

*/

#if !defined(PHEAPPTR_H_4YRJD3YF384Y7RUIHFUFGNXBVBG8423566GJF785673EYGFER7TJG)
#define PHEAPPTR_H_4YRJD3YF384Y7RUIHFUFGNXBVBG8423566GJF785673EYGFER7TJG

#include "pheap.h"
#include "polycast.h"

// the pref_ptr is a descriptor of an object on the persistent heap.
// it has the same(-ish!) semantics as the ref_ptr.

template <class T> class pref_handle;

// pref_ptr

template <class T>
class pref_ptr
{
   public:
      typedef T  element_type;
      typedef T* value_type;

      pref_ptr();

      template <class U> 
      pref_ptr(U* p);

      pref_ptr(pheap::PHeapObject* Obj);

      explicit pref_ptr(pheap::id_type ID);

      pref_ptr(const pref_ptr<T>& r); // copy ctor cannot be a template

      template <class U> pref_ptr(const pref_ptr<U>& r);

      // tagged constructor for dynamic_cast'ing
      struct CastTag {};
      template <class U> pref_ptr(const pref_ptr<U>& r, CastTag);

      ~pref_ptr();

      pref_ptr<T>& operator=(const pref_ptr<T>& r); // assignment operator cannot be a template?
      template <class U> pref_ptr<T>& operator=(const pref_ptr<U>& r);

      pref_ptr<T>& operator=(pheap::PHeapObject* Obj);

      template <class U>
      pref_ptr<T>& operator=(U* Ptr);

      T& operator*() const { DEBUG_PRECONDITION(Ptr != NULL); return *Ptr; }
      T* operator->() const { DEBUG_PRECONDITION(Ptr != NULL); return Ptr; }

      bool operator==(const pref_ptr<T>& r) const { return Ptr == r.Ptr; }
      bool operator!=(const pref_ptr<T>& r) const { return Ptr != r.Ptr; }

      T* get() const { return Ptr; }

      void release(); // releases this.  same as doing *this = pref_ptr();

      pheap::id_type object_id() const { return Handle ? Handle->ID() : 0; }

      pheap::PHeapObject* get_handle() const { return Handle; }

      pref_ptr(T* Ptr_, pheap::PHeapObject* Handle_);

      operator bool() const { return Ptr != NULL; }

      bool is_null() const { return Ptr == NULL; }

   private:
      pheap::PHeapObject* Handle;
      T* Ptr;

   template <class U> friend class pref_ptr;
   template <class U> friend class pref_handle;
};

template <class T, class U>
struct poly_cast_helper<pref_ptr<T>, pref_ptr<U> >
{
   static pref_ptr<T> apply(pref_ptr<U> const& val)
   {
      return pref_ptr<T>(val, pref_ptr<T>::CastTag());
   }
};

template <class T>
inline pref_ptr<T> make_pref_ptr(T* Ptr)
{
   return pref_ptr<T>(Ptr);
}

template <class T>
class pref_handle;

template <int Format, class T>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& stream, pref_handle<T> const& Object);

template <int Format, class T>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& stream, pref_handle<T>& Object);

template <class T>
class pref_handle
{
   public:
      typedef T  element_type;
      typedef T* value_type;

      pref_handle();

      template <class U>
      pref_handle(U* Ptr);

      pref_handle(const pref_handle<T>& r);

      explicit pref_handle(pheap::id_type ID);

      template <class U> pref_handle(const pref_handle<U>& r);
      template <class U> pref_handle(const pref_ptr<U>& r);

      ~pref_handle();

      pref_handle<T>& operator=(const pref_handle<T>& r);

      template <class U> pref_handle<T>& operator=(const pref_handle<U>& r);

      bool operator==(const pref_handle<T>& r) const { return Handle == r.Handle; }
      bool operator!=(const pref_handle<T>& r) const { return Handle != r.Handle; }

      pref_ptr<T> lock() const;
      pref_ptr<T> load() const { return this->lock(); }

      operator bool() const { return Handle != NULL; }

      void release(); // releases this.  same as doing *this = pref_handle();

      pheap::PHeapObject* get_handle() const { return Handle; }

      pheap::id_type object_id() const { return Handle ? Handle->ID() : 0; }

      pref_handle(pheap::PHeapObject* Obj) : Handle(Obj) { }

   private:
      pheap::PHeapObject* Handle;

      template <int Format, typename U>
      friend PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& stream, pref_handle<U> const& Object);
      template <int Format, typename U>
      friend PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& stream, pref_handle<U>& Object);

      template <class U> friend class pref_handle;
};
      
// persistent stream output of pref_handle<T> gives identical results to
// writing a pref_ptr<T> of the same object.
// Similarly, stream input of a pref_handle is interchangeable with pref_ptr.


template <int Format, class T>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& stream, pref_ptr<T> const& Object);

template <int Format, class T>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& stream, pref_ptr<T>& Object);

namespace pheap
{

template <typename T>
void ShutdownPersistent(pref_ptr<T>& MainObject)
{
   PHeapObject* Obj = MainObject.get_handle();
   Obj->AddReference();
   MainObject.release();
   ShutdownPersistent(Obj);
}

template <typename T>
void ShutdownPersistent(pref_handle<T>& MainObject)
{
   PHeapObject* Obj = MainObject.get_handle();
   Obj->AddReference();
   MainObject.release();
   ShutdownPersistent(Obj);
}

template <typename T>
void ExportHeap(PHeapFileSystem::BlockFileSystem* FS, pref_handle<T>& MainObject)
{
   ExportHeap(FS, MainObject.get_handle());
}

template <typename T>
void ExportHeap(PHeapFileSystem::BlockFileSystem* FS, pref_ptr<T>& MainObject)
{
   ExportHeap(FS, MainObject.get_handle());
}

template <typename T>
void ExportHeap(std::string const& FileName, pref_handle<T>& MainObject, 
		int NumFiles = 1, size_t PageSize = 0)
{
   ExportHeap(FileName, MainObject.get_handle(), NumFiles, PageSize);
}

template <typename T>
void ExportHeap(std::string const& FileName, pref_ptr<T>& MainObject, 
		int NumFiles = 1, size_t PageSize = 0)
{
   ExportHeap(FileName, MainObject.get_handle(), NumFiles, PageSize);
}

} // namespace pheap

#include "prefptr.cc"

#endif
