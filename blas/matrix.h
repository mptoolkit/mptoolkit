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
#include "common/types.h"

namespace blas
{

// some future work
class PermutationMatrix;

//
// Memory-based matrix type.  Column-major format for compatability with BLAS.
// Moveable, non-copyable, non-resizable (except by moving).
// Default Tag type is cpu_tag.
//

// forward declaration to construct a dense matrix from a sparse matrix
template <typename T>
class SparseMatrix;

template <typename T, typename Tag>
class Matrix : public NormalMatrix<T, Matrix<T, Tag>, Tag>
{
   public:
      using value_type         = T;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;
      using reference          = typename buffer_type::reference;
      using const_reference    = typename buffer_type::const_reference;

      using NormalMatrix<T, Matrix<T, Tag>, Tag>::operator();

      Matrix() : Rows(0), Cols(0), LeadingDimension(0) {}

      Matrix(Matrix const&) = delete;

      Matrix(Matrix&& Other) = default;

      Matrix(int Rows_, int Cols_, arena Arena_);

      Matrix(int Rows_, int Cols_) : Matrix(Rows_, Cols_, tag_type::template default_arena<T>()) {}

      Matrix(int Rows_, int Cols_, T const& Fill, arena Arena_);

      Matrix(int Rows_, int Cols_, T const& Fill)
         : Matrix(Rows_, Cols_, Fill, tag_type::template default_arena<T>()) {}

      Matrix(SparseMatrix<value_type> const& Other, arena Arena_);

      explicit Matrix(SparseMatrix<value_type> const& Other)
	 : Matrix(Other, tag_type::template default_arena<T>()) {}

      // construction via expression template
      template <typename U>
      Matrix(MatrixRef<T, U, tag_type> const& E, arena Arena_)
         : Matrix(E.rows(), E.cols(), Arena_)
      {
         assign(*this, E.as_derived());
      }

      template <typename U>
      Matrix(MatrixRef<T, U, tag_type> const& E)
         : Matrix(E.rows(), E.cols())
      {
         assign(*this, E.as_derived());
      }

      // construction from intializer list
      template <typename U>
      Matrix(std::initializer_list<std::initializer_list<U>> x, arena Arena_);

      template <typename U>
      Matrix(std::initializer_list<std::initializer_list<U>> x)
         : Matrix(x, tag_type::template default_arena<T>()) {}

      ~Matrix() noexcept = default;

      Matrix& operator=(Matrix const& Other)
      {
	 assign(*this, Other);
	 return *this;
      }

      Matrix& operator=(Matrix&& Other) = default;

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename U>
      Matrix& operator=(MatrixRef<T, U, tag_type> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Matrix& operator+=(MatrixRef<T, U, tag_type> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      Matrix& operator-=(MatrixRef<T, U, tag_type> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      int rows() const { return Rows; }
      int cols() const { return Cols; }

      int leading_dimension() const { return LeadingDimension; }

      static constexpr char trans() { return 'N'; }

      vector_view<T, tag_type>
      row(int r)
      {
         return vector_view<T, tag_type>(Cols, LeadingDimension, Buf.ptr(r));
      }

      const_vector_view<T, tag_type>
      row(int r) const
      {
         return const_vector_view<T, tag_type>(Cols, LeadingDimension, Buf.ptr(r));
      }

      normal_vector_view<T, tag_type>
      column(int c)
      {
         return normal_vector_view<T, tag_type>(Rows, Buf.ptr(LeadingDimension*c));
      }

      const_normal_vector_view<T, tag_type>
      column(int c) const
      {
         return const_normal_vector_view<T, tag_type>(Rows, Buf.ptr(LeadingDimension*c));
      }

      vector_view<T, tag_type>
      diagonal()
      {
         return vector_view<T, tag_type>(std::min(Rows,Cols), LeadingDimension+1, Buf.ptr());
      }

      const_vector_view<T, tag_type>
      diagonal() const
      {
         return const_vector_view<T, tag_type>(std::min(Rows,Cols), LeadingDimension+1, Buf.ptr());
      }

      // sets all elements to zero
      void clear()
      {
         blas::clear(*this);
      }

      buffer_type& buffer() { return Buf; }
      buffer_type const& buffer() const { return Buf; }

      storage_type storage() { return Buf.ptr(); }
      const_storage_type storage() const { return Buf.cptr(); }

      reference operator()(int r, int c)
      {
         DEBUG_RANGE_CHECK(r, 0, Rows);
         DEBUG_RANGE_CHECK(c, 0, Cols);
         return Buf[c*LeadingDimension+r];
      }

      const_reference operator()(int r, int c) const
      {
         DEBUG_RANGE_CHECK(r, 0, Rows);
         DEBUG_RANGE_CHECK(c, 0, Cols);
         return Buf[c*LeadingDimension+r];
      }

      static Matrix make_identity(int Size);

   private:
      int Rows;
      int Cols;
      int LeadingDimension;
      buffer_type Buf;
};

template <typename T, typename Tag>
inline
Matrix<T, Tag>
Matrix<T, Tag>::make_identity(int Size)
{
   Matrix<T, Tag> Result(Size, Size, blas::number_traits<T>::zero());
   fill(Result.diagonal(), blas::number_traits<T>::identity());
   return Result;
}

// copy

template <typename T, typename Tag_Type>
inline
Matrix<T, Tag_Type>
copy(Matrix<T, Tag_Type> const& x, blas::arena const& A)
{
   Matrix<T> Result(x.rows(), x.cols(), A);
   Result = x;
   return Result;
}

template <typename T, typename Tag_Type>
inline
Matrix<T, Tag_Type>
copy(Matrix<T, Tag_Type> const& x)
{
   Matrix<T, Tag_Type> Result(x.rows(), x.cols());
   Result = x;
   return Result;
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

template <typename T, typename Tag>
struct ScalarTypes<blas::Matrix<T, Tag>> : ScalarTypes<T> {};

#include "matrix.icc"

#endif
