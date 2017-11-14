// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix.h
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
// A simple dense matrix class designed for scalar types.
//

#if !defined(MPTOOLKIT_BLAS_MATRIX_H)
#define MPTOOLKIT_BLAS_MATRIX_H

#include "common/trace.h"
#include "arena.h"
#include "matrixref.h"
#include "vector.h"
#include "vector_view.h"
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
class DiagonalMatrix : public DiagonalBlasMatrix<T, DiagonalMatrix<T>, cpu_tag>
{
   public:
      using value_type     = T;
      using iterator       = T*;
      using const_iterator = T const*;
      using pointer        = T*;
      using reference      = T&;

      DiagonalMatrix() = delete;

      DiagonalMatrix(DiagonalMatrix const&) = delete;

      DiagonalMatrix(DiagonalMatrix&& Other) : Size(Other.Size),
			       Arena(std::move(Other.Arena)),
			       Data(Other.Data) { Other.Data = nullptr; }

      DiagonalMatrix(int Rows_, int Cols_, arena Arena_);

      DiagonalMatrix(int Rows_, int Cols_) : DiagonalMatrix(Rows_, Cols_, get_malloc_arena()) {}

      DiagonalMatrix(int Rows_, int Cols_, T const& Fill, arena Arena_);

      DiagonalMatrix(int Rows_, int Cols_, T const& Fill) : DiagonalMatrix(Rows_, Cols_, Fill, get_malloc_arena()) {}

      // construction via expression template
      template <typename U>
      DiagonalMatrix(DiagonalMatrixRef<T, U, cpu_tag> const& E, arena Arena_);

      template <typename U>
      DiagonalMatrix(DiagonalMatrixRef<T, U, cpu_tag> const& E) : DiagonalMatrix(E, get_malloc_arena()) {}

      // construction from intializer list
      template <typename U>
      DiagonalMatrix(std::initializer_list<U> x, arena Arena_);

      template <typename U>
      DiagonalMatrix(std::initializer_list<U> x) : DiagonalMatrix(x, get_malloc_arena()) {}

      ~DiagonalMatrix() { if (Data) Arena.free(Data, Size); }

      DiagonalMatrix& operator=(DiagonalMatrix const& Other)
      {
	 assign(*this, Other);
	 return *this;
      }

      DiagonalMatrix& operator=(DiagonalMatrix&& Other)
      {
	 Size = Other.Size;
	 Arena = std::move(Other.Arena); Data = Other.Data; Other.Data = nullptr;
	 return *this;
      }

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename U>
      DiagonalMatrix& operator=(DiagonalMatrixRef<T, U, cpu_tag> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator+=(DiagonalMatrixRef<T, U, cpu_tag> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator-=(DiagonalMatrixRef<T, U, cpu_tag> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      NormalVectorView<T>
      diagonal()
      {
         return NormalVectorView<T>(Size, Data);
      }

      ConstNormalVectorView<T>
      diagonal() const
      {
         return ConstNormalVectorView<T>(Size, Data);
      }

      int rows() const { return Size; }

      int cols() const { return Size; }

      constexpr int stride() const { return 1; }

      T* storage() { return Data; }
      T const* storage() const { return Data; }

      T& operator()(int i, int ii)
      {
         CHECK_EQUAL(i, ii);
         DEBUG_RANGE_CHECK(i, 0, Size);
         return Data[i];
      }

      T const& operator()(int i, int ii) const
      {
         DEBUG_RANGE_CHECK(i, 0, Size);
         CHECK_EQUAL(i, ii);
         return Data[i];
      }

   private:
      arena Arena;
      int Size;
      T* Data;
};


} // namespace blas

#include "diagonalmatrix.icc"

#endif
