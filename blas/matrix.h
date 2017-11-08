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

// It is important that we include the BLAS and LAPACK headers before
// matrixref.h here so that ADL can find the relevant functions at the
// point of instantiation
#include "common/trace.h"
#include "arena.h"
#include "matrix-blas.h"
#include "matrix-lapack.h"
#include "matrixref.h"
#include "vector.h"
#include "vector_view.h"
#include <list>
#include <mutex>
#include <iostream>
#include <iomanip>
#include "common/formatting.h"
#include "common/randutil.h"

namespace blas
{

// some future work
class PermutationMatrix;

//
// Memory-based matrix type.  Column-major format for compatability with BLAS.
// Moveable, non-copyable, non-resizable (except by moving).
//

template <typename T>
using VectorView = vector_view<T, cpu_tag>;

template <typename T>
using ConstVectorView = const_vector_view<T, cpu_tag>;

template <typename T>
using NormalVectorView = normal_vector_view<T, cpu_tag>;

template <typename T>
using ConstNormalVectorView = const_normal_vector_view<T, cpu_tag>;

template <typename T>
class Matrix : public NormalMatrix<T, Matrix<T>, cpu_tag>
{
   public:
      typedef T value_type;

      Matrix() = delete;

      Matrix(Matrix const&) = delete;

      Matrix(Matrix&& Other) : Rows(Other.Rows), Cols(Other.Cols),
                               LeadingDimension(Other.LeadingDimension),
			       Arena(std::move(Other.Arena)),
			       Data(Other.Data) { Other.Data = nullptr; }

      Matrix(int Rows_, int Cols_, arena Arena_);

      Matrix(int Rows_, int Cols_) : Matrix(Rows_, Cols_, get_malloc_arena()) {}

      Matrix(int Rows_, int Cols_, T const& Fill, arena Arena_);

      Matrix(int Rows_, int Cols_, T const& Fill) : Matrix(Rows_, Cols_, Fill, get_malloc_arena()) {}

      // construction via expression template
      template <typename U>
      Matrix(MatrixRef<T, Matrix<T>, U> const& E, arena Arena_);

      template <typename U>
      Matrix(MatrixRef<T, Matrix<T>, U> const& E) : Matrix(E, get_malloc_arena()) {}

      // construction from intializer list
      template <typename U>
      Matrix(std::initializer_list<std::initializer_list<U>> x, arena Arena_);

      template <typename U>
      Matrix(std::initializer_list<std::initializer_list<U>> x) : Matrix(x, get_malloc_arena()) {}

      template <typename U>
      Matrix(MatrixRef<T, U, cpu_tag> const& E)
         : Matrix(E.rows(), E.cols())
      {
         assign(*this, E.as_derived());
      }

      ~Matrix() { if (Data) Arena.free(Data, LeadingDimension*Cols); }

      Matrix& operator=(Matrix const& Other)
      {
	 assign(*this, Other);
	 return *this;
      }

      Matrix& operator=(Matrix&& Other)
      {
	 Rows = Other.Rows; Cols = Other.Cols; LeadingDimension = Other.LeadingDimension;
	 Arena = std::move(Other.Arena); Data = Other.Data; Other.Data = nullptr;
	 return *this;
      }

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename U>
      Matrix& operator=(MatrixRef<T, U, cpu_tag> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Matrix& operator+=(MatrixRef<T, U, cpu_tag> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Matrix& operator-=(MatrixRef<T, U, cpu_tag> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      int rows() const { return Rows; }
      int cols() const { return Cols; }

      int leading_dimension() const { return LeadingDimension; }

      constexpr char trans() const { return 'N'; }

      VectorView<T>
      row(int r)
      {
         return VectorView<T>(Cols, LeadingDimension, Data + r);
      }

      ConstVectorView<T>
      row(int r) const
      {
         return ConstVectorView<T>(Cols, LeadingDimension, Data + r);
      }

      NormalVectorView<T>
      column(int c)
      {
         return VectorView<T>(Rows, Data + LeadingDimension*c);
      }

      ConstNormalVectorView<T>
      column(int c) const
      {
         return ConstVectorView<T>(Rows, Data + LeadingDimension*c);
      }

      VectorView<T>
      diagonal()
      {
         return VectorView<T>(std::min(Rows,Cols), LeadingDimension+1, Data);
      }

      ConstVectorView<T>
      diagonal() const
      {
         return ConstVectorView<T>(std::min(Rows,Cols), LeadingDimension+1, Data);
      }

      // sets all elements to zero
      void clear()
      {
         std::memset(Data, 0, LeadingDimension*Cols*sizeof(T));
      }

      T* storage() { return Data; }
      T const* storage() const { return Data; }

      T& operator()(int r, int c)
      {
         DEBUG_RANGE_CHECK(r, 0, Rows);
         DEBUG_RANGE_CHECK(c, 0, Cols);
         return Data[c*LeadingDimension+r];
      }

      T const& operator()(int r, int c) const
      {
         DEBUG_RANGE_CHECK(r, 0, Rows);
         DEBUG_RANGE_CHECK(c, 0, Cols);
         return Data[c*LeadingDimension+r];
      }

   private:
      arena Arena;
      int Rows;
      int Cols;
      int LeadingDimension;
      T* Data;
};

// copy

template <typename T>
inline
Matrix<T>
copy(Matrix<T> const& x, blas::arena const& A)
{
   Matrix<T> Result(x.rows(), x.cols(), A);
   Result = x;
   return Result;
}

template <typename T>
inline
Matrix<T>
copy(Matrix<T> const& x)
{
   Matrix<T> Result(x.rows(), x.cols());
   Result = x;
   return Result;
}

template <typename T>
inline
Matrix<T>
copy(Matrix<T>&& x)
{
   return std::move(x);
}

// random_matrix

template <typename T>
Matrix<T>
random_matrix(int Rows, int Cols)
{
   Matrix<T> Result(Rows, Cols);
   for (int r = 0; r < Rows; ++r)
   {
      for (int c = 0; c < Cols; ++c)
      {
         Result(r,c) = randutil::rand<T>();
      }
   }
   return Result;
}

// I/O

template <typename T, typename U>
std::ostream&
operator<<(std::ostream& out, NormalMatrix<T, U, cpu_tag> const& x);

} // namespace blas

#include "matrix.icc"

#endif
