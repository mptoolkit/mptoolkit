// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/vector.h
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
// A simple dense vector class designed for scalar types.
//

#if !defined(MPTOOLKIT_BLAS_VECTOR_H)
#define MPTOOLKIT_BLAS_VECTOR_H

#include "common/trace.h"
#include "arena.h"
#include "vectorref.h"
#include "vector_view.h"
#include <list>
#include <mutex>
#include <iostream>
#include <iomanip>
#include "common/formatting.h"

namespace blas
{

template <typename T>
struct cpu_buffer
{
   arena Arena;
   T*    Ptr;
   int   Size;

   using storage_type       = T*;
   using const_storage_type = T const*;

   using reference       = T&;
   using const_reference = T const&;

   cpu_buffer() = delete;

   cpu_buffer(arena const& a, T* p, int s) : Arena(a), Ptr(p), Size(s) {}
   cpu_buffer(cpu_buffer&& other) : Arena(std::move(other.Arena)), Ptr(other.Ptr), Size(other.Size)
   {
      other.Ptr = nullptr;
   }

   cpu_buffer(cpu_buffer&) = delete;

   ~cpu_buffer()
   {
      if (Ptr)
         Arena.free(Ptr, Size);
      Ptr = nullptr;
   }

   cpu_buffer& operator=(cpu_buffer&) = delete;
   cpu_buffer& operator=(cpu_buffer&&) = delete;

   T* ptr() { return Ptr; }
   T* ptr() const { return Ptr; }
   T const* cptr() const { return Ptr; }

   T* ptr(int offset) { return Ptr + offset; }
   T* ptr(int offset) const { return Ptr + offset; }
   T const* cptr(int offset) const { return Ptr + offset; }

   T& operator[](int offset) { return Ptr[offset]; }
   T const& operator[](int offset) const { return Ptr[offset]; }

   static cpu_buffer allocate(int Size, arena A)
   {
      return cpu_buffer(A, A.allocate_type<T>(Size), Size);
   }
};

struct cpu_tag
{
   template <typename T>
   using buffer_type = cpu_buffer<T>;

   template <typename T>
   using storage_type = T*;

   template <typename T>
   using const_storage_type = T const*;

   template <typename T>
   using async_ref = T&;

   template <typename T>
   static
   arena default_arena() { return get_malloc_arena(); }

   template <typename T>
   static int select_leading_dimension(int ld)
   {
      return ld;
   }
};

//
// Memory-based vector type.
// Moveable, non-copyable, non-resizable (except by moving).
//

template <typename T, typename Tag = cpu_tag>
class Vector : public NormalVector<T, Vector<T>, Tag>
{
   public:
      using value_type     = T;
      using iterator       = T*;
      using const_iterator = T const*;
      using pointer        = T*;
      using reference      = T&;

      Vector() = delete;

      Vector(Vector const&) = delete;

      Vector(Vector&& Other) : Size(Other.Size),
			       Arena(std::move(Other.Arena)),
			       Data(Other.Data) { Other.Data = nullptr; }

      Vector(int Size_, arena Arena_);

      Vector(int Size_) : Vector(Size_, get_malloc_arena()) {}

      Vector(int Size_, T const& Fill, arena Arena_);

      Vector(int Size_, T const& Fill) : Vector(Size_, Fill, get_malloc_arena()) {}

      // construction via expression template
      template <typename U>
      Vector(VectorRef<T, Vector<T>, U> const& E, arena Arena_);

      template <typename U>
      Vector(VectorRef<T, Vector<T>, U> const& E) : Vector(E, get_malloc_arena()) {}

      // construction from intializer list
      template <typename U>
      Vector(std::initializer_list<U> x, arena Arena_);

      template <typename U>
      Vector(std::initializer_list<U> x) : Vector(x, get_malloc_arena()) {}

      ~Vector() { if (Data) Arena.free(Data, Size); }

      Vector& operator=(Vector const& Other)
      {
	 assign(*this, Other);
	 return *this;
      }

      Vector& operator=(Vector&& Other)
      {
	 Size = Other.Size;
	 Arena = std::move(Other.Arena); Data = Other.Data; Other.Data = nullptr;
	 return *this;
      }

      iterator begin() { return Data; }
      iterator end() { return Data + Size; }

      const_iterator begin() const { return Data; }
      const_iterator end() const { return Data + Size; }

      const_iterator cbegin() const { return Data; }
      const_iterator cend() const { return Data + Size; }

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename U>
      Vector& operator=(VectorRef<T, U, cpu_tag> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Vector& operator+=(VectorRef<T, U, cpu_tag> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Vector& operator-=(VectorRef<T, U, cpu_tag> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      int size() const { return Size; }

      constexpr int stride() const { return 1; }

      T* storage() & { return Data; }
      T const* storage() const& { return Data; }

      T& operator[](int i)
      {
         DEBUG_RANGE_CHECK(i, 0, Size);
         return Data[i];
      }

      T const& operator[](int i) const
      {
         DEBUG_RANGE_CHECK(i, 0, Size);
         return Data[i];
      }

   private:
      arena Arena;
      int Size;
      T* Data;
};

template <typename T>
using VectorView = vector_view<T, cpu_tag>;

template <typename T>
using ConstVectorView = const_vector_view<T, cpu_tag>;

template <typename T, typename U>
std::ostream&
operator<<(std::ostream& out, BlasVector<T, U, cpu_tag> const& x)
{
   out << '[' << x.size() << "]\n";
   bool first = true;
   T const* Ptr = x.storage();
   for (int i = 0; i < x.size(); ++i)
   {
      if (!first)
         out << '\n';
      write_format(out, *Ptr);
      Ptr += x.stride();
   }
   return out;
}

// copy

template <typename T, typename U>
Vector<T>
copy(blas::BlasVector<T, Vector<T>, U> const& x, blas::arena const& A)
{
   Vector<T> Result(x.size(), A);
   assign(Result, x.derived());
   return Result;
}

template <typename T, typename U>
Vector<T>
copy(blas::BlasVector<T, Vector<T>, U> const& x)
{
   Vector<T> Result(x.size());
   assign(Result, x.derived());
   return Result;
}

} // namespace blas

#include "vector.icc"

#endif
