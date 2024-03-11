// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/datablock_old.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

size_t const CacheLineSize = 0;  // This is a hook to make the data block aligned on a cache-line.  not yet implemented.

namespace Private
{

// This is a bit nasty; to use a single memory allocation we need to play some tricks
// to get the alignment and offsets of the members.  This relies on using offsetof()
// on a type that isn't strictly a POD.  But in practice the compiler would have to
// go out of its way to break this.
template <typename T>
struct AlignmentHelper
{
   AtomicRefCount RefCount;
   T*             Extent;
   T              Data[1];

};

template <typename T>
struct BlockAlignTraits
{
   static int const ActualRefCountOffset = offsetof(AlignmentHelper<T>, RefCount);
   static int const ActualExtentOffset = offsetof(AlignmentHelper<T>, Extent);
   static int const ActualDataOffset = offsetof(AlignmentHelper<T>, Data);

   // These are relative to the Data field
   static int const RelativeRefCountOffset = ActualRefCountOffset - ActualDataOffset;
   static int const RelativeExtentOffset = ActualExtentOffset - ActualDataOffset;

   // The size we need to allocate is BasicSize + sizeof(T) * number_of_elements
   static int const BasicSize = ActualDataOffset;

   // helper functions to do the nasty offset-shifting stuff
   static T* allocate(size_t Size);

   static void deallocate();

   static AtomicRefCount& GetRefCount(T* Data);
   static T* GetExtent(T* Data);
};

template <typename T>
inline
T* BlockAlignTraits<T>::allocate(size_t Size)
{
   unsigned char* Base = static_cast<unsigned char*>(malloc(BasicSize + sizeof(T) * Size));
   new (Base+ActualRefCountOffset) AtomicRefCount(1);
   new (Base+Actual

} // namespace Private

template <class T>
class DataBlock
{
   public:

      DataBlock() : Data(new T[1]), RefCount(new AtomicRefCount(1)) {}

      explicit DataBlock(size_t Size);

      DataBlock(DataBlock<T> const& Other) : Data(Other.get()), RefCount(Other.RefCount)
      { ++*RefCount; }

      DataBlock& operator=(DataBlock<T> const& Other)
      { ++*Other.RefCount; if (--*RefCount == 0) { delete[] Data; delete RefCount; } Data = Other.Data; RefCount = Other.RefCount;
      return *this; }

      ~DataBlock() { if (--*RefCount == 0) { delete[] Data; delete RefCount; } }

      // returns true if the reference count is > 1, false if reference count == 1.
      bool is_shared() const { return RefCount->value() > 1; }

      T* get() const { return Data; }

   private:
      T* Data;
      AtomicRefCount* RefCount;

   friend class DataBlock<T const>;
};

template <class T>
DataBlock<T>::DataBlock(size_t Size)
   : Data(new T[std::max(size_t(1), Size + CacheLineSize)]), RefCount(new AtomicRefCount(1))
{
}

template <class T>
class DataBlock<T const>
{
   public:
      DataBlock() : Data(new T[1]), RefCount(new AtomicRefCount(1)) {}

      // is this ctor completely useless?
      explicit DataBlock(size_t Size);

      DataBlock(DataBlock<T const> const& Other) : Data(Other.get()), RefCount(Other.RefCount)
      { ++*RefCount; }

      DataBlock(DataBlock<T> const& Other) : Data(Other.get()), RefCount(Other.RefCount)
      { ++*RefCount; }

      DataBlock& operator=(DataBlock<T const> const& Other)
      { ++*Other.RefCount; if (--*RefCount == 0) { delete[] Data; delete RefCount; } Data = Other.Data; RefCount = Other.RefCount;
      return *this; }

      DataBlock& operator=(DataBlock<T> const& Other)
      { ++*Other.RefCount; if (--*RefCount == 0) { delete[] Data; delete RefCount; } Data = Other.Data; RefCount = Other.RefCount;
      return *this; }

      ~DataBlock() { if (--*RefCount == 0) { delete[] Data; delete RefCount; } }

      // returns true if the reference count is > 1, false if reference count == 1.
      bool is_shared() const { return RefCount->value() > 1; }

      T const* get() const { return Data; }

   private:
      T const* Data;
      AtomicRefCount*     RefCount;
};

template <class T>
DataBlock<T const>::DataBlock(size_t Size)
   : Data(new T[std::max(size_t(1), Size + CacheLineSize)]), RefCount(new AtomicRefCount(1))
{
}

#endif
