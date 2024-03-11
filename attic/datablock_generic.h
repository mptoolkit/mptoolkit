// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/datablock_generic.h
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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
/* -*- C++ -*- $Id$
  A fast reference counted block of memory.

  also defines range, which is a template to represent an iterator range, and
  fast_copy, which is similar to std::copy, but only uses preincrement.

  2002-07-01: Fixed a bug in DataBlock::DataBlock(size_t) where allocating size zero would fail.
  It no longer fails, although maybe it should instead be a precondition violation?
*/

#if !defined(DATABLOCK_H_DSHJ34278943UYUJ34OIUR489UWEIO4)
#define DATABLOCK_H_DSHJ34278943UYUJ34OIUR489UWEIO4

#include <stdlib.h>
#include <algorithm>
#include "dataops.h"
#include "common/atomic.h"
#include "common/construct_array.h"

size_t const CacheLineSize = 0;  // This is a hook to make the data block
                                 // aligned on a cache-line.  not yet implemented.

namespace Private
{

template <typename Header>
struct Overhead
{
   AtomicRefCount RefCount;
   size_t         Size;
   Header         Head;

   Overhead(int Count, size_t Size_, Header const& Head_ = Header()) : RefCount(Count), Size(Size_), Head(Head_) {}
};

template <>
struct Overhead<void>
{
   AtomicRefCount RefCount;
   size_t         Size;

   Overhead(int Count, size_t Size_) : RefCount(Count), Size(Size_) {}
};

template <typename T, typename Header>
struct DataBlockAllocator
{
   typedef Private::Overhead<Header> OverheadType;
   static int const ArrayAlign = boost::alignment_of<T>::value;
   static int const ArrayOffset = ((sizeof(OverheadType) - 1) / ArrayAlign + 1) * ArrayAlign;

   static T* allocate(size_t Size);
   static T* copy(T const* Data);
   static void deallocate(T const* Data);

   static AtomicRefCount& GetRefCount(T const* Data);
   static Header& GetHeader(T const* Data);
   static size_t GetSize(T const* Data);
};

template <typename T, typename Header>
inline
T* DataBlockAllocator<T, Header>::allocate(size_t Size)
{
   unsigned char* Base;
   try
   {
      Base = new unsigned char[ArrayOffset + sizeof(T) * Size];
      T* Data = reinterpret_cast<T*>(Base + ArrayOffset);

      new (Base) Overhead<Header>(1, Size);
      ext::construct_array(Data, Size);
      return Data;
   }
   catch (...)
   {
      if (Base) delete[] Base;
      throw;
   }
}

template <typename T, typename Header>
inline
void DataBlockAllocator<T, Header>::deallocate(T const* Data)
{
   unsigned char const* Base = reinterpret_cast<unsigned char const*>(Data) - ArrayOffset;
   size_t Size = reinterpret_cast<Overhead<Header> >(Base)->Size;

   ext::destroy_array(Data, Size);
   reinterpret_cast<Overhead<Header> >(Base)->~Overhead<Header>();
   delete[] Base;
}

template <typename T>
inline
T* DataBlockAllocator<T, Header>::copy(T const* Data)
{
   size_t Size = DataBlockAllocator<T, Header>::GetSize(Data);

   unsigned char* NewBase;
   try
   {
      NewBase = new unsigned char[ArrayOffset + sizeof(T) * Size];
      T* NewData = reinterpret_cast<T*>(NewBase + ArrayOffset);

      new (NewBase) Overhead<Header>(1, Size);
      ext::construct_array(NewData, Size);
      std::uninitialized_copy(Data, Data+Size, NewData);
      return Data;
   }
   catch (...)
   {
      if (Base) delete[] Base;
      throw;
   }
}

template <typename T, typename Header>
inline
AtomicRefCount& DataBlockAllocator<T, Header>::GetRefCount(T const* Data)
{
   unsigned char* Base = reinterpret_cast<unsigned char*>(const_cast<T*>(Data)) - ArrayOffset;
   return reinterpret_cast<Overhead<Header>*>(Base)->RefCount;
}

template <typename T, typename Header>
inline
size_t DataBlockAllocator<T, Header>::GetSize(T* Data)
{
   unsigned char const* Base = reinterpret_cast<unsigned char const*>(const_cast<T*>(Data)) - ArrayOffset;
   return reinterpret_cast<Overhead<Header> const*>(Base)->Size;
}

template <typename T, typename Header>
inline
Header& DataBlockAllocator<T, Header>::GetHeader(T const* Data)
{
   unsigned char* Base = reinterpret_cast<unsigned char*>(const_cast<T*>(Data)) - ArrayOffset;
   return reinterpret_cast<Overhead<Header>*>(Base)->Head;
}

} // namespace Private

template <typename T, typename Header = void>
class DataBlock
{
   private:
      typedef Private::DataBlockAllocator<T, Header> Allocator;

   public:

      DataBlock() : Data(Allocator::allocate(1)) {}

      explicit DataBlock(size_t Size) : Data(Allocator::allocate(Size)) {}

      DataBlock(DataBlock<T> const& Other) : Data(Other.Data)
      { ++this->ref_count(); }

      DataBlock& operator=(DataBlock<T> const& Other)
      { ++Other->ref_count(); if (--this->ref_count() == 0) { Allocator::deallocate(Data); }
      Data = Other.Data;
      return *this; }

      ~DataBlock() { if (--this->ref_count() == 0) { Allocator::deallocate(Data); } }

      // returns true if the reference count is > 1, false if reference count == 1.
      bool is_shared() const { return this->ref_count().value() > 1; }

      // returns a deep copy of the block
      DataBlock<T, Header> deep_copy() const { return DataBlock<T>(Allocator::copy(Data)); }

      void cow() const { if (this->is_shared()) *this = this->deep_copy(); }

      size_t size() const { return Allocator::GetSize(Data); }

      T* get() const { return Data; }

      Header& header() const { return Allocator::GetHeader(Data); }

   private:
      explicit DataBlock(T* D) : Data(D_) {}
      AtomicRefCount& ref_count() const { return Allocator::GetRefCount(Data); }

      char* Data;

   friend class DataBlock<T const, Header>;
};


template <typename T, typename Header>
class DataBlock<T const, Header>
{
   private:
      typedef Private::DataBlockAllocator<T, Header> Allocator;

   public:

      DataBlock() : Data(Allocator::allocate(1)) {}

      explicit DataBlock(size_t Size) : Data(Allocator::allocate(Size)) {}

      DataBlock(DataBlock<T, Header> const& Other) : Data(Other.Data)
      { ++this->ref_count(); }

      DataBlock(DataBlock<T const, Header> const& Other) : Data(Other.Data)
      { ++this->ref_count(); }

      DataBlock& operator=(DataBlock<T, Header> const& Other)
      { ++Other->ref_count(); if (--this->ref_count() == 0) { Allocator::deallocate(Data); }
      Data = Other.Data;
      return *this; }

      DataBlock& operator=(DataBlock<T const, Header> const& Other)
      { ++Other->ref_count(); if (--this->ref_count() == 0) { Allocator::deallocate(Data); }
      Data = Other.Data;
      return *this; }

      ~DataBlock() { if (--this->ref_count() == 0) { Allocator::deallocate(Data); } }

      // returns true if the reference count is > 1, false if reference count == 1.
      bool is_shared() const { return this->ref_count().value() > 1; }

      // returns a deep copy of the block
      DataBlock<T const, Header> deep_copy() const { return DataBlock<T>(Allocator::copy(Data)); }

      size_t size() const { return Allocator::GetSize(Data); }

      T const* get() const { return Data; }

      Header const& header() const { return Allocator::GetHeader(Data); }

   private:
      explicit DataBlock(T const* D) : Data(D_) {}
      AtomicRefCount& ref_count() const { return Allocator::GetRefCount(Data); }

      T const* Data;
};

#endif
