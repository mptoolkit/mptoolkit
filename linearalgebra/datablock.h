// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/datablock.h
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
/* -*- C++ -*- $Id$
  A fast reference counted block of memory.

  2002-07-01: Fixed a bug in DataBlock::DataBlock(size_t) where allocating size zero would fail.
  It no longer fails, although maybe it should instead be a precondition violation?

  2004-04-24: Optimized to use a single memory allocation; the array and the reference count are
              constructed into a single block of uninitialized memory.  This also allows for
              an arbitrary header to be included.  This defaults to void (ie. no extra header).
              The Size is now included in the DataBlock.  This allows DataBlock to implement
              its own deep_copy() and cow() functions for copy-on-write semantics.
              datablockreference.h builds reference semantics on top, via a double indirection scheme.
*/

#if !defined(DATABLOCK_H_DSHJ34278943UYUJ34OIUR489UWEIO4)
#define DATABLOCK_H_DSHJ34278943UYUJ34OIUR489UWEIO4

#include <stdlib.h>
#include <algorithm>
#include "common/atomic.h"
#include "common/construct_array.h"
#include "common/trace.h"
#include <memory>

#if defined(DATABLOCK_TRACE_DETAILED)
#define TRACE_DATABLOCK(Msg) TRACE(Msg)
#else
#define TRACE_DATABLOCK(Msg) DUMMY_TRACE(Msg)
#endif

size_t const CacheLineSize = 0;  // This is a hook to make the data block
                                 // aligned on a cache-line.  not yet implemented.

struct NoHeader {};  // dummy struct for the case where we have no extra header

namespace Private
{

// struct Overhead encapsulates the extra stuff that we construct into the
// same block of memory as the array
template <typename Header>
struct Overhead
{
   AtomicRefCount RefCount;
   size_t         Size;
   Header         Head;

   Header const& get_header() const { return Head; }
   Header& get_header() { return Head; }

   Overhead(int Count, size_t Size_, Header const& Head_)
     : RefCount(Count), Size(Size_), Head(Head_) {}
};

template <typename Dummy = void>
struct NoHeaderDummy
{
  static NoHeader Head;
};

template <typename Dummy>
NoHeader NoHeaderDummy<Dummy>::Head;

template <>
struct Overhead<NoHeader>
{
   AtomicRefCount RefCount;
   size_t         Size;
  //   static NoHeader Head;

   NoHeader& get_header() const { return NoHeaderDummy<>::Head; }

   Overhead(int Count, size_t Size_, NoHeader Head_ = NoHeader()) : RefCount(Count), Size(Size_) {}
};

// DataBlockAllocator is a traits class that does all of the fiddly work of
// manipulating the overhead + array
template <typename T, typename Header>
struct DataBlockAllocator
{
   typedef Private::Overhead<Header> OverheadType;
   // *** ATLAS seems to require alignment of 8 for arrays.  As a work-around,
   // we multiply the alignment by 2.
   static int const ArrayAlign = boost::alignment_of<T>::value * 2;
   static int const ArrayOffset = ((sizeof(OverheadType) - 1) / ArrayAlign + 1) * ArrayAlign;

   static T* create_empty(Header const& Head = Header());
   static T* allocate(size_t Size, Header const& Head = Header());
   static T* allocate_fill(size_t Size, T const& Fill, Header const& Head = Header());
   static T* copy(T const* Data);
   static void deallocate(T const* Data);

   static AtomicRefCount& GetRefCount(T const* Data);
   static Header& GetHeader(T const* Data);
   static size_t GetSize(T const* Data);
};

template <typename T, typename Header>
inline
T* DataBlockAllocator<T, Header>::allocate(size_t Size, Header const& Head)
{
   unsigned char* Base = static_cast<unsigned char*>(::operator new(ArrayOffset + sizeof(T) * Size));
   T* Data = reinterpret_cast<T*>(Base + ArrayOffset);
   try // lame attempt at exception safety
   {
      new (Base) Overhead<Header>(1, Size, Head);
      ext::construct_array(Data, Size);
      return Data;
   }
   catch (...)
   {
      if (Base) ::operator delete(Base);
      throw;
   }
}

template <typename T, typename Header>
inline
T* DataBlockAllocator<T, Header>::create_empty(Header const& Head)
{
   unsigned char* Base = static_cast<unsigned char*>(::operator new(ArrayOffset));
   T* Data = reinterpret_cast<T*>(Base + ArrayOffset);
   try // lame attempt at exception safety
   {
      new (Base) Overhead<Header>(1, 0, Head);
      return Data;
   }
   catch (...)
   {
      if (Base) ::operator delete(Base);
      throw;
   }
}

template <typename T, typename Header>
inline
T* DataBlockAllocator<T, Header>::allocate_fill(size_t Size, T const& Fill, Header const& Head)
{
   unsigned char* Base = static_cast<unsigned char*>(::operator new(ArrayOffset + sizeof(T) * Size));
   T* Data = reinterpret_cast<T*>(Base + ArrayOffset);
   try // lame attempt at exception safety
   {
      new (Base) Overhead<Header>(1, Size, Head);
      ext::construct_array(Data, Size, Fill);
      return Data;
   }
   catch (...)
   {
      if (Base) ::operator delete(Base);
      throw;
   }
}

template <typename T, typename Header>
inline
void DataBlockAllocator<T, Header>::deallocate(T const* Data)
{
   unsigned char const* Base = reinterpret_cast<unsigned char const*>(Data) - ArrayOffset;
   size_t Size = reinterpret_cast<Overhead<Header> const*>(Base)->Size;

   ext::destroy_array(Data, Size);
   reinterpret_cast<Overhead<Header> const*>(Base)->~Overhead<Header>();
   ::operator delete(const_cast<void*>(static_cast<void const*>(Base)));
}

template <typename T, typename Header>
//inline
T* DataBlockAllocator<T, Header>::copy(T const* Data)
{
   size_t Size = DataBlockAllocator<T, Header>::GetSize(Data);

   unsigned char* NewBase = static_cast<unsigned char*>
     (::operator new(ArrayOffset + sizeof(T) * Size));
   T* NewData = reinterpret_cast<T*>(NewBase + ArrayOffset);
   //   try
   //   {
      new (NewBase) Overhead<Header>(1, Size, DataBlockAllocator<T, Header>::GetHeader(Data));
      //      ext::construct_array(NewData, Size);
      std::uninitialized_copy(Data, Data+Size, NewData);
      return NewData;
      //   }
      //   catch (...)
      //   {
      //      ::operator delete(NewBase);
      //      throw;
      //   }
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
size_t DataBlockAllocator<T, Header>::GetSize(T const* Data)
{
   unsigned char const* Base = reinterpret_cast<unsigned char const*>(const_cast<T*>(Data)) - ArrayOffset;
   return reinterpret_cast<Overhead<Header> const*>(Base)->Size;
}

template <typename T, typename Header>
inline
Header& DataBlockAllocator<T, Header>::GetHeader(T const* Data)
{
   unsigned char* Base = reinterpret_cast<unsigned char*>(const_cast<T*>(Data)) - ArrayOffset;
   return reinterpret_cast<Overhead<Header>*>(Base)->get_header();
}

} // namespace Private

template <typename T, typename Header = NoHeader>
class DataBlock
{
   private:
      typedef Private::DataBlockAllocator<T, Header> Allocator;

   public:

      DataBlock() : Data(Allocator::create_empty()) { TRACE_DATABLOCK("DataBlock: default constructing"); }

      explicit DataBlock(size_t Size, Header const& Head = Header())
        : Data(Allocator::allocate(Size, Head))
      { TRACE_DATABLOCK("DataBlock: constructing")(Size) << "Data = " << (void*) Data; }

      explicit DataBlock(Header const& Head)
        : Data(Allocator::create_empty(Head))
      { TRACE_DATABLOCK("DataBlock: constructing"); }

      DataBlock(size_t Size, T const& Fill, Header const& Head = Header())
        : Data(Allocator::allocate_fill(Size, Fill, Head))
      { TRACE_DATABLOCK("DataBlock: constructing")(Size) << "Data = " << (void*) Data; }

      DataBlock(DataBlock const& Other) : Data(Other.Data)
      { ++this->ref_count();
        TRACE_DATABLOCK("DataBlock: copy constructing")(this->ref_count().value())
        << "Other.Data = " << (void*) Other.Data; }

      DataBlock& operator=(DataBlock const& Other)
      { TRACE_DATABLOCK("DataBlock: assignment") << "Other.Data = " << (void*) Other.Data;
      ++Other.ref_count(); this->do_sub_reference();
      Data = Other.Data;
      return *this; }

      ~DataBlock() { this->do_sub_reference(); }

      AtomicRefCount& ref_count() const { return Allocator::GetRefCount(Data); }

      // returns true if the reference count is > 1, false if reference count == 1.
      bool is_shared() const { return this->ref_count().value() > 1; }

      // returns a deep copy of the block
      DataBlock<T, Header> deep_copy() const { return DataBlock<T, Header>(Allocator::copy(Data)); }

      void cow();

      size_t size() const { return Allocator::GetSize(Data); }

      T* get() const { return Data; }

      Header& header() const { return Allocator::GetHeader(Data); }

   private:
      void do_sub_reference();
      explicit DataBlock(T* D_) : Data(D_) {}

      T* Data;

   friend class DataBlock<T const, Header>;
};


template <typename T, typename Header>
inline
void DataBlock<T, Header>::cow()
{
   TRACE_DATABLOCK("DataBlock::cow()")(this->is_shared()) << "Data = " << (void*) Data;
   if (this->is_shared())
   {
      T* NewData = Allocator::copy(Data);
      this->do_sub_reference();
      Data = NewData;
      TRACE_DATABLOCK("DataBlock::cow()") << "New Data = " << (void*) Data;
   }
}


template <typename T, typename Header>
inline
void DataBlock<T, Header>::do_sub_reference()
{
   if (--this->ref_count() == 0)
   {
      TRACE_DATABLOCK("DataBlock::do_sub_reference: deallocating") << "Data = " << (void*) Data;
      Allocator::deallocate(Data);
   }
   else
   {
      TRACE_DATABLOCK("DataBlock::do_sub_reference: reference count remains positive.");
   }
}

// DataBlock<T const> probably isn't useful anymore

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

      Header const* header() const { return Allocator::GetHeader(Data); }

   private:
      explicit DataBlock(T const* D_) : Data(D_) {}
      AtomicRefCount& ref_count() const { return Allocator::GetRefCount(Data); }

      T const* Data;
};

#endif
