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

#if !defined(MPTOOLKIT_BLAS_DIAGONALMATRIX_H)
#define MPTOOLKIT_BLAS_DIAGONALMATRIX_H

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

template <typename T, typename Tag = cpu_tag>
class DiagonalMatrix : public DiagonalBlasMatrix<T, DiagonalMatrix<T, Tag>, Tag>
{
   public:
      using value_type         = T;
      using tag_type           = Tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;
      using reference          = typename buffer_type::reference;
      using const_reference    = typename buffer_type::const_reference;

      DiagonalMatrix() = delete;

      DiagonalMatrix(DiagonalMatrix const&) = delete;

      DiagonalMatrix(DiagonalMatrix&& Other) = default;

      DiagonalMatrix(int Rows_, int Cols_, arena Arena_);

      DiagonalMatrix(int Rows_, int Cols_) : DiagonalMatrix(Rows_, Cols_, tag_type::template default_arena<T>()) {}

      DiagonalMatrix(int Rows_, int Cols_, T const& Fill, arena Arena_);

      DiagonalMatrix(int Rows_, int Cols_, T const& Fill)
         : DiagonalMatrix(Rows_, Cols_, Fill, tag_type::template default_arena<T>()) {}

      DiagonalMatrix(int Rows_) : DiagonalMatrix(Rows_, Rows_) {}

      // construction via expression template
      template <typename U>
      DiagonalMatrix(DiagonalMatrixRef<T, U, tag_type> const& E, arena Arena_);

      template <typename U>
      DiagonalMatrix(DiagonalMatrixRef<T, U, tag_type> const& E) : DiagonalMatrix(E, tag_type::template default_arena<T>()) {}

      // construction from intializer list
      template <typename U>
      DiagonalMatrix(std::initializer_list<U> x, arena Arena_);

      template <typename U>
      DiagonalMatrix(std::initializer_list<U> x) : DiagonalMatrix(x, tag_type::template default_arena<T>()) {}

      ~DiagonalMatrix() noexcept = default;

      DiagonalMatrix& operator=(DiagonalMatrix const& Other)
      {
	 assign(*this, Other);
	 return *this;
      }

      DiagonalMatrix& operator=(DiagonalMatrix&& Other) = default;

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename U>
      DiagonalMatrix& operator=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator+=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator-=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      normal_vector_view<T, tag_type>
      diagonal()
      {
         return normal_vector_view<T, tag_type>(Size, Buf.ptr());
      }

      const_normal_vector_view<T, tag_type>
      diagonal() const
      {
         return const_normal_vector_view<T, tag_type>(Size, Buf.cptr());
      }

      int rows() const { return Size; }

      int cols() const { return Size; }

      constexpr int stride() const { return 1; }

      buffer_type& buffer() { return Buf; }
      buffer_type const& buffer() const { return Buf; }

      storage_type storage() { return Buf.ptr(); }
      const_storage_type storage() const { return Buf.cptr(); }

      reference operator()(int i, int ii)
      {
         CHECK_EQUAL(i, ii);
         DEBUG_RANGE_CHECK(i, 0, Size);
         return Buf[i];
      }

      const_reference operator()(int i, int ii) const
      {
         CHECK_EQUAL(i, ii);
         DEBUG_RANGE_CHECK(i, 0, Size);
         return Buf[i];
      }

      static DiagonalMatrix make_identity(int Size);

   private:
      int Size;
      buffer_type Buf;
};

template <typename T>
class DiagonalMatrix<T, cpu_tag> : public DiagonalBlasMatrix<T, DiagonalMatrix<T, cpu_tag>, cpu_tag>
{
   public:
      using value_type         = T;
      using iterator           = T*;
      using const_iterator     = T const*;
      using tag_type           = cpu_tag;
      using buffer_type        = typename tag_type::template buffer_type<T>;
      using storage_type       = typename buffer_type::storage_type;
      using const_storage_type = typename buffer_type::const_storage_type;
      using reference          = typename buffer_type::reference;
      using const_reference    = typename buffer_type::const_reference;
      using row_type           = T;

      DiagonalMatrix() = delete;

      DiagonalMatrix(DiagonalMatrix const&) = delete;

      DiagonalMatrix(DiagonalMatrix&& Other) = default;

      DiagonalMatrix(int Rows_, int Cols_, arena Arena_);

      DiagonalMatrix(int Rows_, int Cols_) : DiagonalMatrix(Rows_, Cols_, tag_type::template default_arena<T>()) {}

      DiagonalMatrix(int Rows_, int Cols_, T const& Fill, arena Arena_);

      DiagonalMatrix(int Rows_, int Cols_, T const& Fill)
         : DiagonalMatrix(Rows_, Cols_, Fill, tag_type::template default_arena<T>()) {}

      explicit DiagonalMatrix(int Rows_) : DiagonalMatrix(Rows_, Rows_) {}

      // construction via expression template
      template <typename U>
      DiagonalMatrix(DiagonalMatrixRef<T, U, tag_type> const& E, arena Arena_);

      template <typename U>
      DiagonalMatrix(DiagonalMatrixRef<T, U, tag_type> const& E) : DiagonalMatrix(E, tag_type::template default_arena<T>()) {}

      // construction from intializer list
      template <typename U>
      DiagonalMatrix(std::initializer_list<U> x, arena Arena_);

      template <typename U>
      DiagonalMatrix(std::initializer_list<U> x) : DiagonalMatrix(x, tag_type::template default_arena<T>()) {}

      ~DiagonalMatrix() noexcept = default;

      DiagonalMatrix& operator=(DiagonalMatrix const& Other)
      {
	 assign(*this, Other);
	 return *this;
      }

      DiagonalMatrix& operator=(DiagonalMatrix&& Other) = default;

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename U>
      DiagonalMatrix& operator=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator+=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator-=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      normal_vector_view<T, tag_type>
      diagonal()
      {
         return normal_vector_view<T, tag_type>(Size, Buf.ptr());
      }

      const_normal_vector_view<T, tag_type>
      diagonal() const
      {
         return const_normal_vector_view<T, tag_type>(Size, Buf.cptr());
      }

      int rows() const { return Size; }

      int cols() const { return Size; }

      constexpr int stride() const { return 1; }

      buffer_type& buffer() { return Buf; }
      buffer_type const& buffer() const { return Buf; }

      storage_type storage() { return Buf.ptr(); }
      const_storage_type storage() const { return Buf.cptr(); }

      reference operator()(int i, int ii)
      {
         CHECK_EQUAL(i, ii);
         DEBUG_RANGE_CHECK(i, 0, Size);
         return Buf[i];
      }

      const_reference operator()(int i, int ii) const
      {
         CHECK_EQUAL(i, ii);
         DEBUG_RANGE_CHECK(i, 0, Size);
         return Buf[i];
      }

      template <typename U>
      void insert(int r, int c, U const& value)
      {
	 DEBUG_CHECK_EQUAL(r,c);
         DEBUG_RANGE_CHECK_OPEN(r, 0, Size);
	 Buf[r] = value;
      }

      template <typename U>
      void insert(int r, int c, U&& value)
      {
	 DEBUG_CHECK_EQUAL(r,c);
         DEBUG_RANGE_CHECK_OPEN(r, 0, Size);
	 Buf[r] = std::move(value);
      }

      template <typename U>
      void add(int r, int c, U const& value)
      {
	 DEBUG_CHECK_EQUAL(r,c);
         DEBUG_RANGE_CHECK_OPEN(r, 0, Size);
	 Buf[r] += value;
      }

      template <typename U>
      void add(int r, int c, U&& value)
      {
	 DEBUG_CHECK_EQUAL(r,c);
         DEBUG_RANGE_CHECK_OPEN(r, 0, Size);
	 Buf[r] += std::move(value);
      }

      template <typename U>
      void subtract(int r, int c, U const& value)
      {
	 DEBUG_CHECK_EQUAL(r,c);
         DEBUG_RANGE_CHECK_OPEN(r, 0, Size);
	 Buf[r] -= value;
      }

      template <typename U>
      void subtract(int r, int c, U&& value)
      {
	 DEBUG_CHECK_EQUAL(r,c);
         DEBUG_RANGE_CHECK_OPEN(r, 0, Size);
	 Buf[r] -= std::move(value);
      }

      static DiagonalMatrix make_identity(int Size);

      iterator begin() { return Buf.ptr(); }
      iterator end() { return Buf.ptr() + Size; }

      const_iterator begin() const { return Buf.cptr(); }
      const_iterator end() const { return Buf.cptr() + Size; }

      const_iterator cbegin() const { return Buf.cptr(); }
      const_iterator cend() const { return Buf.cptr() + Size; }

   private:
      int Size;
      buffer_type Buf;
};

template <typename T, typename Tag>
inline
DiagonalMatrix<T, Tag>
DiagonalMatrix<T, Tag>::make_identity(int Size)
{
   DiagonalMatrix<T, Tag> Result(Size, Size, blas::number_traits<T>::identity());
   return Result;
}

} // namespace blas

#include "diagonalmatrix.icc"

#endif
