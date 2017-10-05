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
#include <list>
#include <mutex>
#include <iostream>
#include <iomanip>
#include "common/formatting.h"

namespace blas
{

// some future work
class PermutationMatrix;

//
// Memory-based matrix type.  Column-major format for compatability with BLAS.
// Moveable, non-copyable, non-resizable (except by moving).
//

template <typename T>
class Matrix;

template <typename T>
struct blas_traits<Matrix<T>>
{
   using storage_type       = T*;
   using const_storage_type = T const*;
};

template <typename T>
class Matrix : public BlasMatrix<T, Matrix<T>>
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

      int leading_dimension() const { return LeadingDimension; }

      constexpr char trans() const { return 'N'; }

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

template <typename T>
void
wite_format(std::ostream& out, T x)
{
   out << std::setw(6) << x;
}

void
write_format(std::ostream& out, double x)
{
   out << std::setw(10) << x;
}

void
write_format(std::ostream& out, std::complex<double> x)
{
   out << format_complex(x);
}

template <typename T>
std::ostream&
operator<<(std::ostream& out, Matrix<T> const& x)
{
   out << '[' << x.rows() << ", " << x.cols() << "]\n";
   for (int r = 0; r < x.rows(); ++r)
   {
      bool first = true;
      for (int c = 0; c < x.cols(); ++c)
      {
         if (!first)
            out << ", ";
         write_format(out, x(r,c));
      }
      out << '\n';
   }
   return out;
}

} // namespace blas

#include "matrix.icc"

#endif
