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
#include <list>
#include <mutex>
#include <iostream>
#include <iomanip>
#include "common/formatting.h"

namespace blas
{

//
// Memory-based vector type.
// Moveable, non-copyable, non-resizable (except by moving).
//

template <typename T>
class Vector;

template <typename T>
struct blas_traits<Vector<T>>
{
   using storage_type       = T*;
   using const_storage_type = T const*;
};

template <typename T>
class Vector : public BlasVector<T, Vector<T>>
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
      Vector& operator=(VectorRef<T, Vector<T>, U> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Vector& operator+=(VectorRef<T, Vector<T>, U> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Vector& operator-=(VectorRef<T, Vector<T>, U> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      int size() const { return Size; }

      constexpr int stride() const { return 1; }

      T* storage() { return Data; }
      T const* storage() const { return Data; }

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
std::ostream&
operator<<(std::ostream& out, Vector<T> const& x)
{
   out << '[' << x.size() << "]\n";
   bool first = true;
   for (auto const& a : x)
   {
      if (!first)
         out << '\n';
      write_format(out, a);
      first = false;
   }
   return out;
}

} // namespace blas

#include "vector.icc"

#endif
