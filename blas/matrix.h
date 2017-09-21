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

#include <list>
#include <mutex>
#include <iostream>
#include "arena.h"
#include "matrixref.h"

namespace blas
{

// some future work
class PermutationMatrix;

//
// Memory-based matrix type.  Column-major format for compatability with BLAS.
// Moveable, non-copyable, non-resizable (except by moving).
//

template <typename T>
class Matrix : public BlasMatrix<T, Matrix<T>>
{
   public:
      typedef T value_type;

      Matrix() = delete;

      Matrix(Matrix const&) = delete;

      Matrix(Matrix&& Other) : Rows(Other.Rows), Cols(Other.Cols), LeadingDimension(Other.LeadingDimension),
			       Arena(std::move(Other.Arena)),
			       Data(Other.Data) { Other.Data = nullptr; }

      Matrix(int Rows_, int Cols_, arena Arena_ = MallocArena);

      Matrix(int Rows_, int Cols_, T const& Fill, arena Arena_ = MallocArena);

      // construction via expression template
      template <typename U>
      Matrix(MatrixRef<T, Matrix<T>, U> const& E, arena Arena_ = MallocArena);

      ~Matrix() { if (Data) Arena.free(Data, LeadingDimension*Cols); }

      Matrix& operator=(Matrix const& Other)
      {
	 assign(*this, Other);
	 return *this;
      }

      Matrix& operator=(Matrix&& Other)
      {
	 Rows = Outer.Rows; Cols = Other.Cols; LeadingDimension = Other.LeadingDimension;
	 Arena = std::move(Other.Arens); Data = Other.Data; Other.Data = nullptr;
	 return *this;
      }

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename V>
      Matrix& operator=(MatrixRef<T, Matrix<T>, U> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Matrix& operator+=(MatrixRef<T, Matrix<T>, U> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Matrix& operator-=(MatrixRef<T, Matrix<T>, U> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      int rows() const { return Rows; }
      int cols() const { return Cols; }

      int leading_dim() const { return LeadingDimension; }

      char trans() constexpr { return 'N'; }

      T* data() { return Buf.device_ptr(); }
      T const* data() const { return Buf.device_ptr(); }

   private:
      arena Arena;
      int Rows;
      int Cols;
      int LeadingDimension;
      T* Data;
};

} // namespace blas

#endif
