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

namespace detail
{

struct value_ptr_control_block
{
   value_ptr_control_block() = delete;

   value_ptr_control_block(std::function<void()> d) noexcept
      : Deleter(std::move(d)), UseCount(1)
   {
   }

   template <typename T>
   value_ptr_control_block(T* ptr, std::function<void(T*)> d) noexcept
   : value_ptr_control_block([]{ d(ptr); })
   {
   }

   template <typename T>
   value_ptr_control_block(T* ptr) noexcept
   : value_ptr_control_block([]{ delete ptr; })
   {
   }

   ~value_ptr_control_block() noexcept = default;

   void add_ref() noexcept { ++UseCount; }
   int sub_ref() noexcept { return --UseCount; }

   bool unique() const noexcept { return UseCount.is_unique(); }

   int use_count() const noexcept { return UseCount.value(); }

   void invoke_deleter() const noexcept
   {
      Deleter();
   }

   std::function<void()> Deleter;
   AtomicRefCount UseCount;
   // for the ensure_unique() function we need to track the number
   // of threads that are currently making a copy.  (or take the simple
   // route and always make a copy if use_count > 1, in which case two threads
   // that simultaneously do a copy_on_write will generate two new copies, ultimately
   // destroying the original object.  One of these copies is unnecessary.
};

} // namespace detail
