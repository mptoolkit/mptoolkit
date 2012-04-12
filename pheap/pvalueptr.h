/* -*- C++ -*- $Id$

  Created 2000-1-17
  Ian McCulloch

  A persistent version of the value_ptr class.

  2000-06-06: Cleaned up a bit, pvalue_lock added

*/

#if !defined(PVALUEPOINTER_H_FA768SF687SF6FDSH2U9C898Y7GVRJUIGE5U89)
#define PVALUEPOINTER_H_FA768SF687SF6FDSH2U9C898Y7GVRJUIGE5U89

#include "pheap.h"
#include "polycast.h"

// the pvalue_ptr is a descriptor of an object on the persistent heap, with 
// value semantics.

// pvalue_lock.  The idea of this class is that it avoids multiple calls to pvalue_ptr::mutate()
// when makeing repeated accesses to a pvalue_ptr.  Usage would be something like
//
// pvalue_ptr<T> foo = bar;
// pvalue_ptr<T>::lock foo_lock(foo_lock);
// while (cond) 
//    foo_lock->modify();
//
// this is equivalent to
//
// while (cond)
//    foo.mutate()->modify();
//
// except that the overhead of mutate() checking to see if it needs to do copy-on-write
// at every call is eliminated.

template <class T>
class pvalue_handle;

template <class T>
class pvalue_ptr;

template <class T>
class pvalue_lock
{
   public:
      explicit pvalue_lock(pvalue_ptr<T>&);

      T& operator*()             { return *Ptr; }
      T const& operator*() const { return *Ptr; }

      T* operator->()             { return Ptr; }
      T const* operator->() const { return Ptr; }

   private:
      pvalue_lock(); // not implemented
      explicit pvalue_lock(T* Ptr_) : Ptr(Ptr_) {}

      T* Ptr;
};

template <class T>
class pvalue_ptr
{
   public:
      typedef T              value_type;
      typedef T*             pointer;
      typedef T const*       const_pointer;
      typedef T&             reference;
      typedef T const&       const_reference;
      typedef pvalue_lock<T> lock_type;

      pvalue_ptr();

      pvalue_ptr(const pvalue_ptr<T>& r); // copy ctor cannot be a template
      template <class U> pvalue_ptr(const pvalue_ptr<U>& r);

      ~pvalue_ptr();

      template <class U>
      pvalue_ptr(U* ptr);

      pvalue_ptr(pheap::Private::PHeapObject* obj);

      explicit pvalue_ptr(pheap::id_type ID);

      // tagged constructor for dynamic_cast'ing
      struct CastTag {};
      template <class U> pvalue_ptr(const pvalue_ptr<U>& r, CastTag);

      pvalue_ptr<T>& operator=(const pvalue_ptr<T>& r); // assignment operator cannot be a template?

      template <class U> pvalue_ptr<T>& operator=(const pvalue_ptr<U>& r);

      const T& operator*() const { PRECONDITION(Ptr != NULL); return *Ptr; }
      const T* operator->() const { PRECONDITION(Ptr != NULL); return Ptr; }

      T* mutate() { PRECONDITION(Ptr != NULL); cow(); return Ptr; }

      lock_type lock() { return lock_type(*this); }

      bool operator==(const pvalue_ptr<T>& r) const { return Ptr == r.Ptr; }
      bool operator!=(const pvalue_ptr<T>& r) const { return Ptr != r.Ptr; }

      pheap::id_type object_id() const { return Handle ? Handle->ID() : 0; }

      bool is_null() const { return Ptr == NULL; }

      const T* get() const { return Ptr; }

      void release(); // releases this.  same as doing *this = pvalue_ptr();

      pheap::Private::PHeapObject* get_handle() const { return Handle; }

      //      operator bool() const { return Ptr != NULL; }

   private:
      pvalue_ptr(T* Ptr_, pheap::Private::PHeapObject* Handle_) : Handle(Handle_), Ptr(Ptr_) {}

      pheap::Private::PHeapObject* Handle;
      T* Ptr;

      void cow();

   template <class U> friend class pvalue_ptr;
   template <class U> friend class pvalue_handle;

   // SGI compiler is broken & wont accept explicit function template parameter.
   // This doesnt need to be a template friend.
   template <class U> friend pvalue_ptr<U> make_pvalue_ptr(U* Ptr);

   template <class U> friend PStream::ipstream& operator>>(PStream::ipstream& in, pvalue_ptr<U>& Object);
   template <class U> friend PStream::opstream& operator<<(PStream::opstream& out, pvalue_ptr<U> Object);
};

template <typename T>
bool operator<(pvalue_ptr<T> const& x, pvalue_ptr<T> const& y); // not implemented

template <class T>
class pvalue_handle;

template <class T>
PStream::opstream& operator<<(PStream::opstream& stream, pvalue_handle<T> Object);

template <class T>
PStream::ipstream& operator>>(PStream::ipstream& stream, pvalue_handle<T>& Object);

template <class T>
class pvalue_handle
{
   public:
      typedef T  element_type;
      typedef T* value_type;
      typedef pvalue_ptr<T> ptr_type;

      pvalue_handle();

      template <class U>
      pvalue_handle(U* Ptr);

      pvalue_handle(const pvalue_handle<T>& r);

      explicit pvalue_handle(pheap::id_type ID);

      template <class U> pvalue_handle(const pvalue_handle<U>& r);
      template <class U> pvalue_handle(const pvalue_ptr<U>& r);

      ~pvalue_handle();

      pvalue_handle<T>& operator=(const pvalue_handle<T>& r);

      template <class U> pvalue_handle<T>& operator=(const pvalue_handle<U>& r);

      bool operator==(const pvalue_handle<T>& r) const { return Handle == r.Handle; }
      bool operator!=(const pvalue_handle<T>& r) const { return Handle != r.Handle; }

      pvalue_ptr<T> lock() const;
      pvalue_ptr<T> load() const { return this->lock(); }

      bool is_null() const { return Handle == NULL; }

      //      operator bool() const { return Handle != NULL; }

      void release(); // releases this.  same as doing *this = pvalue_handle();

      pheap::PHeapObject* get_handle() const { return Handle; }

      pheap::id_type object_id() const { return Handle ? Handle->ID() : 0; }

      pvalue_handle(pheap::PHeapObject* Obj) : Handle(Obj) { }

   private:
      pheap::PHeapObject* Handle;

      friend PStream::opstream& operator<< <T>(PStream::opstream& stream, pvalue_handle<T> Object);
      friend PStream::ipstream& operator>> <T>(PStream::ipstream& stream, pvalue_handle<T>& Object);

      template <class U> friend class pvalue_handle;
};
      
template <typename T>
bool operator<(pvalue_handle<T> const& x, pvalue_handle<T> const& y); // not implemented

// persistent stream output of pvalue_handle<T> gives identical results to
// writing a pvalue_ptr<T> of the same object.
// Similarly, stream input of a pvalue_handle is interchangeable with pvalue_ptr.

template <class T, class U>
struct poly_cast_helper<pvalue_ptr<T>, pvalue_ptr<U> >
{
   static pvalue_ptr<T> apply(pvalue_ptr<U> const& val)
   {
      return pvalue_ptr<T>(val, pvalue_ptr<T>::CastTag());
   }
};

// inlines

template <class T>
inline
pvalue_lock<T>::pvalue_lock(pvalue_ptr<T>& PVal) : Ptr(PVal.mutate())
{
}

template <class T>
inline pvalue_ptr<T> make_pvalue_ptr(T* Ptr)
{
   return pvalue_ptr<T>(Ptr);
}

template <class T>
inline void pvalue_ptr<T>::cow()
{
   PRECONDITION(Handle != NULL);
   Handle = Handle->CopyOnWrite(Ptr);
   Handle->SetDirty();
}

namespace pheap
{

template <typename T>
void ShutdownPersistent(pvalue_ptr<T>& MainObject)
{
   PHeapObject* Obj = MainObject.get_handle();
   Obj->AddReference();
   MainObject.release();
   ShutdownPersistent(Obj);
}
template <typename T>
void ExportHeap(std::string const& FileName, pvalue_ptr<T>& MainObject,
		int NumFiles = 1, size_t PageSize = 0)
{
   ExportHeap(FileName, MainObject.get_handle(), NumFiles, PageSize);
}

} // namespace pheap

#include "pvalueptr.cc"

#endif

