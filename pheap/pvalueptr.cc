// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pvalueptr.cc
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

template <class T>
pvalue_ptr<T>::pvalue_ptr() : Handle(NULL), Ptr(NULL) {}

template <typename T>
pvalue_ptr<T>::pvalue_ptr(pheap::Private::PHeapObject* Obj)
  : Handle(Obj), Ptr(Obj->template Lock<T>())
{
   Handle->SubReference();
}

template <typename T>
pvalue_ptr<T>::pvalue_ptr(pheap::id_type ID)
{
   std::tie(Ptr, Handle) = pheap::GetObject<T>(ID);
}

template <class T>
template <class U>
pvalue_ptr<T>::pvalue_ptr(U* p) : Handle(p ? pheap::Allocate(p) : NULL), Ptr(p)
{
}

template <class T>
pvalue_ptr<T>::pvalue_ptr(const pvalue_ptr<T>& r) : Handle(r.Handle), Ptr(r.Ptr)
{
   if (Handle) Handle->AddLock();
}

template <class T>
template <class U>
pvalue_ptr<T>::pvalue_ptr(const pvalue_ptr<U>& r) : Handle(r.Handle), Ptr(r.Ptr)
{
   if (Handle) Handle->AddLock();
}

template <class T>
template <class U>
inline
pvalue_ptr<T>::pvalue_ptr(const pvalue_ptr<U>& r, CastTag) : Handle(r.Handle),Ptr(dynamic_cast<T*>(r.Ptr))
{
   if (Ptr == NULL) Handle = NULL;
   else Handle->AddLock();
}

template <class T>
pvalue_ptr<T>& pvalue_ptr<T>::operator=(const pvalue_ptr<T>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.Handle->AddLock();
   if (Handle) Handle->SubLock();
   Handle = r.Handle;
   Ptr = r.Ptr;
   return *this;
}

template <class T>
template <class U>
pvalue_ptr<T>& pvalue_ptr<T>::operator=(const pvalue_ptr<U>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.Handle->AddLock();
   if (Handle) Handle->SubLock();
   Handle = r.Handle;
   Ptr = r.Ptr;
   return *this;
}

template <class T>
void pvalue_ptr<T>::release()
{
   if (Handle) Handle->SubLock();
   Handle = NULL;
   Ptr = NULL;
}

template <class T>
pvalue_ptr<T>::~pvalue_ptr()
{
   if (Handle) Handle->SubLock();
}

template <class T>
PStream::opstream& operator<<(PStream::opstream& stream, pvalue_ptr<T> Object)
{
   stream << "$vp";
   stream.put_id(Object.object_id());
   return stream;
}

template <class T>
PStream::ipstream& operator>>(PStream::ipstream& stream, pvalue_ptr<T>& Object)
{
   std::string uid;
   stream >> uid;
   if (!(uid == "$vp")) { PANIC("bad pstream, expecting a uid of $vp")(uid); }
   Object = pvalue_ptr<T>(stream.get_id());
   return stream;
}

// pvalue_handle

template <class T>
inline
pvalue_handle<T>::pvalue_handle() : Handle(NULL)
{
}

template <class T>
template <class U>
inline
pvalue_handle<T>::pvalue_handle(U* Ptr)
{
   if (Ptr)
   {
      T* CvtPtr = static_cast<T*>(Ptr);
      Handle = PHeapFileSystem::Allocate(CvtPtr);
      Handle->AddReference();
      Handle->SubLock();
   }
   else Handle = NULL;
}

template <class T>
pvalue_handle<T>::pvalue_handle(pheap::id_type ID)
  : Handle(pheap::AddReference(ID))
{
}

template <class T>
inline
pvalue_handle<T>::pvalue_handle(const pvalue_handle<T>& r) : Handle(r.Handle)
{
   if (Handle) Handle->AddReference();
}

template <class T>
inline
pvalue_handle<T>::pvalue_handle(const pvalue_ptr<T>& r) : Handle(r.Handle)
{
   if (Handle) Handle->AddReference();
}

#if 0
template <class T>
template <class U>
inline
pvalue_handle<T>::pvalue_handle(const pvalue_handle<U>& r) : Handle(r.Handle)
{
   if (Handle) Handle->AddReference();
}

template <class T>
template <class U>
inline
pvalue_handle<T>::pvalue_handle(const pvalue_ptr<U>& r) : Handle(r.Handle)
{
   if (Handle) Handle->AddReference();
}
#endif

template <class T>
inline
pvalue_handle<T>::~pvalue_handle()
{
   release();
}

template <class T>
inline
pvalue_handle<T>&
pvalue_handle<T>::operator=(const pvalue_handle<T>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.Handle->AddReference();
   if (Handle) Handle->SubReference();
   Handle = r.Handle;
   return *this;
}
#if 0
template <class T>
template <class U>
inline
pvalue_handle<T>&
pvalue_handle<T>::operator=(const pvalue_handle<U>& r)
{
   // it is essential that we increment r's reference count before we do any delete
   // so that assignment of a pointer to the same object doesnt delete the object
   if (r.Handle) r.Handle->AddReference();
   if (Handle) Handle->SubReference();
   Handle = r.Handle;
   return *this;
}
#endif

template <class T>
pvalue_ptr<T> pvalue_handle<T>::lock() const
{
   if (!Handle) return pvalue_ptr<T>();
   T* Ptr = Handle->template Lock<T>();

   // Patched, 2004-12-06, needs verifying.
   // We really don't want to SetDirty here, pvalue_ptr::cow() will do that
   // only if we request a non-const pointer.
   //   Handle->SetDirty();

   return pvalue_ptr<T>(Ptr, Handle);
}

template <class T>
inline
void
pvalue_handle<T>::release()
{
   if (Handle) Handle->SubReference();
   Handle = NULL;
}

template <class T>
PStream::opstream& operator<<(PStream::opstream& stream, pvalue_handle<T> Object)
{
   stream << "$vp";
   stream.put_id(Object.object_id());
   return stream;
}

template <class T>
PStream::ipstream& operator>>(PStream::ipstream& stream, pvalue_handle<T>& Object)
{
   std::string uid;
   stream >> uid;
   if (!(uid == "$vp")) { PANIC("bad stream uid, expecting $pp")(uid); }
   Object = pvalue_handle<T>(stream.get_id());
   return stream;
}
