// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/prefptr.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "pheap.h"
#include "common/bindpair.h"

#include <iostream>
#include <typeinfo>  // for debugging

template <class T>
pref_ptr<T>::pref_ptr() : Handle(NULL), Ptr(NULL)
{
}

template <typename T>
pref_ptr<T>::pref_ptr(pheap::PHeapObject* Obj)
  : Handle(Obj), Ptr(Obj->template Lock<T>())
{
   Handle->SetDirty();
   Handle->SubReference();
}

template <typename T>
pref_ptr<T>::pref_ptr(T* Ptr_, pheap::PHeapObject* Obj)
  : Handle(Obj), Ptr(Ptr_)
{
   // the caller should already have done Handle->SetDirty() here
}

template <typename T>
pref_ptr<T>::pref_ptr(pheap::id_type ID)
{
   std::tie(Ptr, Handle) = pheap::GetObject<T>(ID);
   if (Handle) Handle->SetDirty();
}

template <class T>
template <class U>
pref_ptr<T>::pref_ptr(U* p) : Handle(p ? PHeapFileSystem::Allocate(p) : NULL), Ptr(p)
{
   //      TRACE("printing")(Handle);
   //      Handle->DebugPrint(std::cerr);
}

template <class T>
inline
pref_ptr<T>::pref_ptr(pref_ptr<T> const& r) : Handle(r.Handle), Ptr(r.Ptr)
{
   if (Handle) Handle->AddLock();
}

template <class T>
template <class U>
inline
pref_ptr<T>::pref_ptr(pref_ptr<U> const& r) : Handle(r.Handle), Ptr(r.get())
{
   if (Handle) Handle->AddLock();
}

template <class T>
template <class U>
inline
pref_ptr<T>::pref_ptr(const pref_ptr<U>& r, CastTag) : Handle(r.Handle), Ptr(dynamic_cast<T*>(r.get()))
{
   if (Ptr == NULL) Handle = NULL;
   else Handle->AddLock();
}

template <class T>
template <class U>
pref_ptr<T>& pref_ptr<T>::operator=(U* P)
{
   release();

   Handle = P ? PHeapFileSystem::Allocate(P) : NULL;
   Ptr = P;
   return *this;
}

template <class T>
pref_ptr<T>& pref_ptr<T>::operator=(const pref_ptr<T>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.Handle->AddLock();
   if (Handle) Handle->SubLock();
   Handle = r.Handle;
   Ptr = r.Ptr;
   return *this;
}

template <typename T>
pref_ptr<T>& pref_ptr<T>::operator=(pheap::PHeapObject* Obj)
{
   if (Obj)
   {
      T* P = Obj->template Lock<T>();
      Obj->SubReference();
      if (Handle) Handle->SubLock();
      Handle = Obj;
      Ptr = P;
   }
   else
   {
      Handle = NULL;
      Ptr = NULL;
   }
   return *this;
}


template <class T>
template <class U>
pref_ptr<T>& pref_ptr<T>::operator=(const pref_ptr<U>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.get_handle()->AddLock();
   if (Handle) Handle->SubLock();
   Handle = r.Handle;
   Ptr = r.Ptr;
   return *this;
}

template <class T>
void pref_ptr<T>::release()
{
   if (Handle) Handle->SubLock();
   Handle = NULL;
   Ptr = NULL;
}

template <class T>
pref_ptr<T>::~pref_ptr()
{
   if (Handle) Handle->SubLock();
}


// pref_handle

template <class T>
inline
pref_handle<T>::pref_handle() : Handle(NULL)
{
}

template <class T>
template <class U>
inline
pref_handle<T>::pref_handle(U* Ptr)
{
   if (Ptr)
   {
      T* CvtPtr = static_cast<T*>(Ptr);
      Handle = PHeapFileSystem::Allocate(CvtPtr);
      TRACE("printing")(Handle);
      Handle->DebugPrint(std::cerr);
      Handle->AddReference();
      Handle->SubLock();
   }
   else Handle = NULL;
}

template <class T>
pref_handle<T>::pref_handle(pheap::id_type ID)
  : Handle(pheap::AddReference(ID))
{
}

template <class T>
inline
pref_handle<T>::pref_handle(const pref_handle<T>& r) : Handle(r.Handle)
{
   if (Handle) Handle->AddReference();
}

template <class T>
template <class U>
inline
pref_handle<T>::pref_handle(const pref_handle<U>& r) : Handle(r.Handle)
{
   if (Handle) Handle->AddReference();
}

template <class T>
template <class U>
inline
pref_handle<T>::pref_handle(const pref_ptr<U>& r) : Handle(r.Handle)
{
      TRACE("printing")(Handle);
      Handle->DebugPrint(std::cerr);
   if (Handle) Handle->AddReference();
      TRACE("printing")(Handle);
      Handle->DebugPrint(std::cerr);
}

template <class T>
inline
pref_handle<T>::~pref_handle()
{
   release();
}

template <class T>
inline
pref_handle<T>&
pref_handle<T>::operator=(const pref_handle<T>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.Handle->AddReference();
   if (Handle) Handle->SubReference();
   Handle = r.Handle;
   return *this;
}

template <class T>
template <class U>
inline
pref_handle<T>&
pref_handle<T>::operator=(const pref_handle<U>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.Handle->AddReference();
   if (Handle) Handle->SubReference();
   Handle = r.Handle;
   return *this;
}

template <class T>
pref_ptr<T> pref_handle<T>::lock() const
{
   TRACE("lock()")(Handle);
   if (!Handle) return pref_ptr<T>();
   T* Ptr = Handle->template Lock<T>();
   Handle->SetDirty();
   return pref_ptr<T>(Ptr, Handle);
}

template <class T>
inline
void
pref_handle<T>::release()
{
   if (Handle) Handle->SubReference();
   Handle = NULL;
}

template <int Format, class T>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& stream, pref_handle<T> const& Object)
{
   stream << "$pp";
   stream.put_id(Object.object_id());
   return stream;
}

template <int Format, class T>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& stream, pref_ptr<T> const& Object)
{
   stream << "$pp";
   stream.put_id(Object.object_id());
   return stream;
}

template <int Format, class T>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& stream, pref_ptr<T>& Object)
{
   std::string uid;
   stream >> uid;
   if (!(uid == "$pp")) PANIC("bad stream uid, expecting $pp")(uid);
   Object = pref_ptr<T>(stream.get_id());
   return stream;
}

template <int Format, class T>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& stream, pref_handle<T>& Object)
{
   std::string uid;
   stream >> uid;
   if (!(uid == "$pp")) PANIC("bad stream uid, expecting $pp")(uid);
   Object = pref_handle<T>(stream.sget_id());
   return stream;
}
