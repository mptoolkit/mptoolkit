// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/atomicrefcount.h
//
// Copyright (C) 2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  value_ptr

  A smart pointer with immutable or copy-on-write value sematics.

  The interface mostly mirrors std::shared_ptr, but adds copy-on-write
  semantics via the unique() and ensure_unique() functions.  Unfortunately
  it isn't possible to implement value_ptr using shared_ptr for one
  annoying reason: the use_count() function doesn't provide atomic
  synchronization so it does not have a well-defined value in multi-threaded
  programs.  Similarly, the shared_ptr<T>::unique() function doesn't have
  a well-defined value, and this function is deprecated in c++17.

  For copy-on-write semantics, the required synchronization is minimal: the
  unique() function needs to use acquire semantics when reading the use_count.
  This is sufficient because in a multithread environment, having two threads
  simultaneously reading (eg to copy the value) and writing to the same
  value_ptr<T> is erroneous anyway, so unique() doesn't have to worry about races
  with other threads incrementing the reference count beyond 1 at the same time
  that we check on the reference count.  But we do need to handle the case where
  multiple value_ptr<T> that share the representation simultaneously write (and force
  a copy-on-write).

  value_ptr doesn't support array types.
*/

#if !defined(MPTOOLKIT_COMMON_VALUE_PTR_H)
#define MPTOOLKIT_COMMON_VALUE_PTR_H

#include "atomicrefcount.h"
#include "detail/value_ptr_control.h"
#include "trace.h"
#include <memory>
#include <ostream>

template <typename T>
class value_ptr
{
   public:
      using element_type = T;

      value_ptr() noexcept;
      value_ptr(std::nullptr_t) noexcept;

      template <typename U>
      value_ptr(U* ptr);

      // Construction with a custom deleter.  d(ptr) must be a valid
      // expression.
      template <typename U, typename Deleter>
      value_ptr(U* ptr, Deleter d);

      // Aliasing constructor: share ownership information with x,
      // but get() always returns ptr.  Typically used where ptr
      // is a sub-object of x.
      template <typename U>
      value_ptr(value_ptr<U> const& x, element_type* ptr) noexcept;

      // copy ctor
      value_ptr(value_ptr const& x) noexcept;

      // move ctor
      template <typename U>
      value_ptr(value_ptr<U> const& x) noexcept;

      value_ptr(value_ptr&& x) noexcept;

      template <typename U>
      value_ptr(value_ptr<U>&& x) noexcept;

      // move-construction from a std::unique_ptr
      template< class U, class Deleter>
      value_ptr(std::unique_ptr<U, Deleter>&& r);

      ~value_ptr();

      // copy-assignment
      value_ptr& operator=(value_ptr const& x) noexcept;

      template <typename U>
      value_ptr& operator=(value_ptr<U> const& x) noexcept;

      // move-assignment
      value_ptr& operator=(value_ptr&& x) noexcept;

      template <typename U>
      value_ptr& operator=(value_ptr<U>&& x) noexcept;

      // move assignment from a std::unique_ptr
      template <typename U, typename Deleter>
      value_ptr& operator=(std::unique_ptr<U, Deleter>&& x) noexcept;

      void reset() noexcept;

      template <typename U>
      void reset(U* ptr);

      template <typename U, typename Deleter>
      void reset(U* ptr, Deleter d);

      void swap(value_ptr& x) noexcept;

      // gets the object, returns a const reference
      T const* get() const noexcept { return Ptr; }

      // returns true if the object pointer is not shared;
      // this means that the object can be modified freely
      bool unique() const noexcept;

      // precondition: unique()
      // returns a modifiable pointer
      T* assert_unique() const noexcept;

      // If the use_count is > 1, then force a copy-on-write,
      // and return a modifiable pointer
      T* ensure_unique();

      T const& operator*() const noexcept ( DEBUG_CHECK(Ptr); return *Ptr; }
      T const* operator->() const noexcept { return Ptr; }

      // we don't provide non-const versions of operator* or operator-> since
      // it is too easy to call these by accident (in principle they should be
      // const member functions anyway), so mutating accesses should go via
      // assert_unique() or ensure_unique().

      // returns the number of different value_ptr instances managing the current
      // object.  This is for debugging purposes only, in a multithread environment
      // the value is not well defined, not even to determine if use_count() == 1.
      // Use unique() to determine properly if this is the only instance managing the
      // object.
      int use_count() const noexcept;

      explicit operator bool() const noexcept { return Ptr != nullptr; }

   private:
      T* Ptr;
      detail::value_ptr_control_block* Control;

};

template<class T, class... Args>
value_ptr<T> make_value_ptr(Args&&... args);

template <class T, class U>
inline
bool operator==(value_ptr<T>& lhs, value_ptr<U> const& rhs) noexcept
{
   return lhs.get() == rhs.get();
}

template <class T, class U>
inline
bool operator!=(value_ptr<T>& lhs, value_ptr<U> const& rhs) noexcept
{
   return lhs.get() != rhs.get();
}

template <class T, class U>
inline
bool operator<(value_ptr<T>& lhs, value_ptr<U> const& rhs) noexcept
{
   return lhs.get() < rhs.get();
}

template <class T, class U, class V>
std::basic_ostream<U, V>&
operator<<(std::basic_ostream<U, V>& os, value_ptr<T> const& ptr)
{
   os << ptr.get();
   return os;
}

#include "value_ptr.icc"

#endif
