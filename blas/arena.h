// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/arena.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

//
// Memory allocation arena.
// This is needed for GPU memory, but might also be useful
// for specialized allocators for main memory.
//

#if !defined(MPTOOLKIT_BLAS_ARENA_H)
#define MPTOOLKIT_BLAS_ARENA_H

#include <memory>

namespace blas
{

// rounding up is a useful function for memory allocations...
inline
std::size_t round_up(std::size_t numToRound, std::size_t multiple)
{
    if (multiple == 0)
        return numToRound;

    std::size_t remainder = numToRound % multiple;
    if (remainder == 0)
        return numToRound;

    return numToRound + multiple - remainder;
}

// base class for allocators

class AllocatorBase
{
   public:
      AllocatorBase() = default;

      virtual void* allocate(std::size_t Size, std::size_t Align) = 0;

      virtual void* allocate(std::size_t Size) = 0;

      virtual void free(void* Ptr, std::size_t Size) = 0;

      virtual ~AllocatorBase() noexcept = 0;
};

// memory arena

class arena
{
   public:
      arena() {}

      explicit arena(std::shared_ptr<AllocatorBase> Alloc_) : Alloc(Alloc_) {}

      explicit arena(AllocatorBase* Alloc_) : Alloc(Alloc_) {}

      ~arena() noexcept = default;

      template <typename T>
      T* allocate_type(std::size_t Size) const
      {
         return static_cast<T*>(this->allocate(Size*sizeof(T), sizeof(T)));
      }

      void* allocate(std::size_t Size, std::size_t Align) const
      { return Alloc->allocate(Size, Align); }

      void* allocate(std::size_t Size) const
      { return Alloc->allocate(Size); }

      void free(void* Ptr, std::size_t Size) const
      {
         if (Ptr)
            Alloc->free(Ptr, Size);
      }

   private:
      std::shared_ptr<AllocatorBase> Alloc;
};

class MallocAllocator : public AllocatorBase
{
   public:
      MallocAllocator() {}

      virtual void* allocate(std::size_t Size, std::size_t Align)
      {
         return std::malloc(Size);
      }

      virtual void* allocate(std::size_t Size)
      {
         return std::malloc(Size);
      }

      virtual void free(void* Ptr, std::size_t Size)
      {
         std::free(Ptr);
      }

      virtual ~MallocAllocator() noexcept {}
};

inline
arena const& get_malloc_arena()
{
   static arena A(new MallocAllocator());
   return A;
}

// A 'default' arena for use in allocations involving objects of type T.
// To use, specialize for some type and add an Arena static member.
template <typename T>
struct default_arena
{
};

} // namespace blas

#endif
