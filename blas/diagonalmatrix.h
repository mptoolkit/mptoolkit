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

template <typename T, typename Tag>
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

      DiagonalMatrix() noexcept 
      : Size(0), Buf{} {}

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
	 assign(this->diagonal(), Other.diagonal());
	 return *this;
      }

      DiagonalMatrix& operator=(DiagonalMatrix&& Other) noexcept = default;

      // assignment of expressions based on the same matrix type -- we don't allow assignment
      // of expression templates of other matrix types (eg gpu_matrix)
      // We could allow assignment of different concrete types, but probably not advisable,
      // eg they would generally block, and we can get the same effect with
      // move construction/assignment and get/set operations.
      template <typename U>
      DiagonalMatrix& operator=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 assign(this.diagonal(), E.as_derived().diagonal());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator+=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 add(this->diagonal(), E.as_derived().diagonal());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator-=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 subtract(this->diagonal(), E.as_derived().diagonal());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator*=(U const& x)
      {
	 scale(this->diagonal(), x);
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

//
// a DiagonalMatrix that exists ono the CPU.  This has direct access to memory,
// and an iterator interface.
//

// proxy type for a an element of a diagonal matrix row.
// This is the 'inner' iterator of a 2-dimensional iteration scheme.

template <typename T>
struct DiagonalMatrixElementRef
{
   int Col;
   T& value;

   DiagonalMatrixElementRef() = delete;

   DiagonalMatrixElementRef(int Col_, T& value_)
	 : Col(Col_), value(value_) {}
   
   int col() const { return Col; }

   operator T&() const { return value; }
};

//iterator to a column of a DiagonalMatrixRow

template <typename T>
class DiagonalMatrixColIterator
{
   public:
      using value_type = DiagonalMatrixElementRef<T>;
      using pointer    = T*;
      using reference  = DiagonalMatrixElementRef<T>;
      using iterator_category = std::random_access_iterator_tag;

      DiagonalMatrixColIterator() noexcept
      : Col(0), Ptr(nullptr) {}

      DiagonalMatrixColIterator(int Col_, T* Ptr_)
	 : Col(Col_), Ptr(Ptr_) {}

      int col() const { return Col; }

      reference operator*() const
      {
	 return reference(Col, *Ptr);
      }

      pointer operator->() const
      {
	 return Ptr;
      }

      reference operator[](int i) const
      {
	 return reference(Col+i, Ptr[i]);
      }

      DiagonalMatrixColIterator& operator++() noexcept
      {
	 ++Col;
	 ++Ptr;
	 return *this;
      }

      DiagonalMatrixColIterator operator++(int) noexcept
      {
	 return DiagonalMatrixColIterator(Col++, Ptr++);
      }

      DiagonalMatrixColIterator& operator--() noexcept
      {
	 --Col;
	 --Ptr;
	 return *this;
      }

      DiagonalMatrixColIterator operator--(int) noexcept
      {
	 return DiagonalMatrixColIterator(Col--, Ptr--);
      }

      bool operator==(DiagonalMatrixColIterator const& x) const
      {
	 return Ptr == x.Ptr;
      }

      bool operator!=(DiagonalMatrixColIterator const& x) const
      {
	 return Ptr != x.Ptr;
      }

   private:
      int Col;
      T* Ptr;
};

template <typename T>
class ConstDiagonalMatrixColIterator
{
   public:
      using value_type = DiagonalMatrixElementRef<T const>;
      using pointer    = T const*;
      using reference  = DiagonalMatrixElementRef<T const>;
      using iterator_category = std::random_access_iterator_tag;

      ConstDiagonalMatrixColIterator() noexcept
      : Col(0), Ptr(nullptr) {}

      ConstDiagonalMatrixColIterator(int Col_, T const* Ptr_)
	 : Col(Col_), Ptr(Ptr_) {}

      int col() const { return Col; }

      reference operator*() const
      {
	 return reference(Col, *Ptr);
      }

      pointer operator->() const
      {
	 return Ptr;
      }

      reference operator[](int i) const
      {
	 return reference(Col+i, Ptr[i]);
      }

      ConstDiagonalMatrixColIterator& operator++() noexcept
      {
	 ++Col;
	 ++Ptr;
	 return *this;
      }

      ConstDiagonalMatrixColIterator operator++(int) noexcept
      {
	 return ConstDiagonalMatrixColIterator(Col++, Ptr++);
      }

      ConstDiagonalMatrixColIterator& operator--() noexcept
      {
	 --Col;
	 --Ptr;
	 return *this;
      }

      ConstDiagonalMatrixColIterator operator--(int) noexcept
      {
	 return ConstDiagonalMatrixColIterator(Col--, Ptr--);
      }

      bool operator==(ConstDiagonalMatrixColIterator const& x) const
      {
	 return Ptr == x.Ptr;
      }

      bool operator!=(ConstDiagonalMatrixColIterator const& x) const
      {
	 return Ptr != x.Ptr;
      }

   private:
      int Col;
      T const* Ptr;
};


// proxy type for a row of a diagonal matrix.  This contains only one element,
// and exists primarily for consistent iterator semantics with a SparseMatrix

template <typename T>
class DiagonalMatrixRowRef
{
   public:
      using iterator       = DiagonalMatrixColIterator<T>;
      using const_iterator = ConstDiagonalMatrixColIterator<T>;

      DiagonalMatrixRowRef(int Row_, T& value_)
	 : Row(Row_), value(value_) {}

      // we can point one beyond the end because the value will always be a member of an
      // array

      iterator begin() noexcept { return iterator(Row, &value); }
      const_iterator begin() const noexcept { return const_iterator(Row, &value); }
      const_iterator cbegin() const noexcept { return const_iterator(Row, &value); }

      iterator end() noexcept { return iterator(Row+1, &value+1); }
      const_iterator end() const noexcept { return const_iterator(Row+1, &value+1); }
      const_iterator cend() const noexcept { return const_iterator(Row+1, &value+1); }

      int row() const { return Row; }

      T& operator()(int c)
      {
	 DEBUG_CHECK_EQUAL(c, Row);
	 return value;
      }

      T const& operator()(int c) const
      {
	 DEBUG_CHECK_EQUAL(c, Row);
	 return value;
      }

      iterator find(int Col)
      {
	 if (Col == Row)
	    return iterator(Row, &value);
	 return iterator(Row, &value+1);
      }

      const_iterator find(int Col) const
      {
	 if (Col == Row)
	    return const_iterator(Row, &value);
	 return const_iterator(Row+1, &value+1);
      }

   private:
      int Row;
      T& value;
      
};

template <typename T>
class DiagonalMatrixRowIterator
{
   public:
      using value_type = DiagonalMatrixRowRef<T>;
      using pointer    = T*;
      using reference  = DiagonalMatrixRowRef<T>;
      using iterator_category = std::random_access_iterator_tag;

      DiagonalMatrixRowIterator() noexcept
      : Row(0), Ptr(nullptr) {}

      DiagonalMatrixRowIterator(int Row_, T* Ptr_)
	 : Row(Row_), Ptr(Ptr_) {}

      int row() const { return Row; }

      reference operator*() const
      {
	 return reference(Row, *Ptr);
      }

      pointer operator->() const
      {
	 return Ptr;
      }

      reference operator[](int i) const
      {
	 return reference(Row+i, Ptr[i]);
      }

      DiagonalMatrixRowIterator& operator++() noexcept
      {
	 ++Row;
	 ++Ptr;
	 return *this;
      }

      DiagonalMatrixRowIterator operator++(int) noexcept
      {
	 return DiagonalMatrixRowIterator(Row++, Ptr++);
      }

      DiagonalMatrixRowIterator& operator--() noexcept
      {
	 --Row;
	 --Ptr;
	 return *this;
      }

      DiagonalMatrixRowIterator operator--(int) noexcept
      {
	 return DiagonalMatrixRowIterator(Row--, Ptr--);
      }

      bool operator==(DiagonalMatrixRowIterator const& x) const
      {
	 return Ptr == x.Ptr;
      }

      bool operator!=(DiagonalMatrixRowIterator const& x) const
      {
	 return Ptr != x.Ptr;
      }

   private:
      int Row;
      T* Ptr;
};

template <typename T>
class ConstDiagonalMatrixRowIterator
{
   public:
      using value_type = DiagonalMatrixRowRef<T const>;
      using pointer    = T const*;
      using reference  = DiagonalMatrixRowRef<T const>;
      using iterator_category = std::random_access_iterator_tag;

      ConstDiagonalMatrixRowIterator() noexcept
      : Row(0), Ptr(nullptr) {}

      ConstDiagonalMatrixRowIterator(int Row_, T const* Ptr_)
	 : Row(Row_), Ptr(Ptr_) {}

      int row() const { return Row; }

      reference operator*() const
      {
	 return reference(Row, *Ptr);
      }

      pointer operator->() const
      {
	 return Ptr;
      }

      reference operator[](int i) const
      {
	 return reference(Row+i, Ptr[i]);
      }

      ConstDiagonalMatrixRowIterator& operator++() noexcept
      {
	 ++Row;
	 ++Ptr;
	 return *this;
      }

      ConstDiagonalMatrixRowIterator operator++(int) noexcept
      {
	 return ConstDiagonalMatrixRowIterator(Row++, Ptr++);
      }

      ConstDiagonalMatrixRowIterator& operator--() noexcept
      {
	 --Row;
	 --Ptr;
	 return *this;
      }

      ConstDiagonalMatrixRowIterator operator--(int) noexcept
      {
	 return ConstDiagonalMatrixRowIterator(Row--, Ptr--);
      }

      bool operator==(ConstDiagonalMatrixRowIterator const& x) const
      {
	 return Ptr == x.Ptr;
      }

      bool operator!=(ConstDiagonalMatrixRowIterator const& x) const
      {
	 return Ptr != x.Ptr;
      }

   private:
      int Row;
      T const* Ptr;
};

template <typename T>
class DiagonalMatrix<T, cpu_tag> : public DiagonalBlasMatrix<T, DiagonalMatrix<T, cpu_tag>, cpu_tag>
{
   public:
      using value_type          = T;
      using tag_type            = cpu_tag;
      using buffer_type         = typename tag_type::template buffer_type<T>;
      using storage_type        = typename buffer_type::storage_type;
      using const_storage_type  = typename buffer_type::const_storage_type;
      using reference           = typename buffer_type::reference;
      using const_reference     = typename buffer_type::const_reference;

      using row_reference       = DiagonalMatrixRowRef<T>;
      using const_row_reference = DiagonalMatrixRowRef<T const>;
      using iterator            = DiagonalMatrixRowIterator<T>;
      using const_iterator      = ConstDiagonalMatrixRowIterator<T>;

      DiagonalMatrix() noexcept
      : Size(0) {}

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
	 assign(this->diagonal(), Other.diagonal());
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
	 assign(this->diagonal(), E.as_derived().diagonal());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator+=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 add(this->diagonal(), E.as_derived().diagonal());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator-=(DiagonalMatrixRef<T, U, tag_type> const& E)
      {
	 subtract(this->diagonal(), E.as_derived().diagonal());
	 return *this;
      }

      template <typename U>
      DiagonalMatrix& operator*=(U const& x)
      {
	 scale(*this, x);
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

      row_reference row(int r)
      {
	 DEBUG_RANGE_CHECK(r, 0, Size);
	 return row_reference(r, Buf[r]);
      }

      const_row_reference row(int r) const
      {
	 DEBUG_RANGE_CHECK(r, 0, Size);
	 return const_row_reference(r, Buf[r]);
      }

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

      iterator begin() noexcept { return iterator(0, Buf.ptr()); }
      iterator end() noexcept { return iterator(Size, Buf.ptr() + Size); }

      const_iterator begin() const noexcept { return const_iterator(0, Buf.ptr()); }
      const_iterator end() const noexcept { return const_iterator(Size, Buf.ptr() + Size); }

      const_iterator cbegin() const noexcept { return const_iterator(0, Buf.ptr()); }
      const_iterator cend() const noexcept { return const_iterator(Size, Buf.ptr() + Size); }

   private:
      int Size;
      buffer_type Buf;
};

template <typename T, typename Tag>
DiagonalMatrix<T, Tag>
copy(DiagonalMatrix<T, Tag> const& x)
{
   DiagonalMatrix<T, Tag> Result(x.rows());
   Result = x;
   return Result;
}

// conj

template <typename T>
inline
void
inplace_conj(DiagonalMatrix<T>& x)
{
   inplace_conj(x.diagonal());
}

template <typename T>
inline
DiagonalMatrix<T>
conj(DiagonalMatrix<T> const& x)
{
   using std::conj;
   DiagonalMatrix<T> Result(x.rows(), x.cols());
   Result.diagonal() = conj(x.diagonal());
   return Result;
}

template <typename T>
inline
DiagonalMatrix<T>
conj(DiagonalMatrix<T>&& x)
{
   inplace_conj(x);
   return std::move(x);
}

template <typename T, typename U, typename Tag>
inline
void
norm_frob_sq(blas::DiagonalBlasMatrix<T, U, Tag> const& x,
	     typename Tag::template async_ref<decltype(norm_frob(std::declval<T>()))>& z)
{
   // TODO: DiagonalBlasMatrix should implement its own diagonal() function
   norm_frob_sq(x.as_derived().diagonal(), z);
}

template <typename T, typename U, typename Tag>
//decltype(norm_frob_sq(std::declval<T>()))
auto norm_frob_sq(blas::DiagonalBlasMatrix<T, U, Tag> const& x)
{
   using blas::norm_frob_sq;
   typename Tag::template async_ref<decltype(norm_frob_sq(std::declval<T>()))> z;
   norm_frob_sq(x, z);
   return get_wait(z);
}

template <typename T, typename U, typename Tag>
auto norm_frob(blas::DiagonalBlasMatrix<T, U, Tag> const& x)
{
   using std::sqrt;
   return sqrt(norm_frob_sq(x));
}

template <typename T, typename Tag>
inline
DiagonalMatrix<T, Tag>
DiagonalMatrix<T, Tag>::make_identity(int Size)
{
   DiagonalMatrix<T, Tag> Result(Size, Size, blas::number_traits<T>::identity());
   return Result;
}

template <typename T, typename U, typename Tag, typename V>
inline
void scale(DiagonalBlasMatrix<T, U, Tag>& M, V const& Factor)
{
   // TODO: DiagonalBlasMatrix should have its own diagonal() function
   scale(M.as_derived().diagonal(), Factor);
}

// numeric_type_of specialization
template <typename T, typename Tag>
struct numeric_type_of<DiagonalMatrix<T, Tag>> : numeric_type_of<T> {};

template <typename T, typename Tag>
std::ostream&
operator<<(std::ostream& out, DiagonalMatrix<T, Tag> const& x);

template <int Format, typename T, typename Tag>
PStream::opstreambuf<Format>&
operator<<(PStream::opstreambuf<Format>& out, DiagonalMatrix<T, Tag> const& M);

} // namespace blas

namespace PStream
{

template <typename T, typename Tag>
struct stream_extract<blas::DiagonalMatrix<T, Tag>>
{
   template <int Format>
   blas::DiagonalMatrix<T,Tag> operator()(PStream::ipstreambuf<Format>& in) const;
};
template <typename T>
struct stream_extract<blas::DiagonalMatrix<T, blas::cpu_tag>>
{
   template <int Format>
   blas::DiagonalMatrix<T,blas::cpu_tag> operator()(PStream::ipstreambuf<Format>& in) const;
};

} // namespace PStream

#include "diagonalmatrix.icc"

#endif
