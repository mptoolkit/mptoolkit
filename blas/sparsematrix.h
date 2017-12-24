// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/sparsematrix.h
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
// A generic sparse matrix class
//

#if !defined(MPTOOLKIT_BLAS_SPARSEMATRIX_H)
#define MPTOOLKIT_BLAS_SPARSEMATRIX_H

#include "common/trace.h"
#include "detail/sparserowiterator.h"
#include "number_traits.h"
#include <vector>
#include <map>

namespace blas
{

// moveable but not copyable sparse vector type, based around std::map.

template <typename T, typename IteratorType, typename ConstIteratorType>
class BasicSparseVector
{
   private:
      using base_iterator       = typename std::map<int, T>::iterator;
      using const_base_iterator = typename std::map<int, T>::const_iterator;

   public:
      using value_type       = T;
      using iterator         = IteratorType;
      using const_iterator   = ConstIteratorType;


      BasicSparseVector() = default;
      BasicSparseVector(BasicSparseVector&& Other) = default;
      BasicSparseVector& operator=(BasicSparseVector&& Other) = default;
      ~BasicSparseVector() = default;

      BasicSparseVector(BasicSparseVector&) = delete;
      BasicSparseVector& operator=(BasicSparseVector&) = delete;

      int nnz() const { return Elements.size(); }

      iterator begin() noexcept { return iterator(Elements.begin()); }
      const_iterator begin() const noexcept { return const_iterator(Elements.begin()); }
      const_iterator cbegin() const noexcept { return const_iterator(Elements.cbegin()); }

      iterator end() noexcept { return iterator(Elements.end()); }
      const_iterator end() const noexcept { return const_iterator(Elements.end()); }
      const_iterator cend() const noexcept { return const_iterator(Elements.cend()); }

      iterator find(int Col)
      {
         return iterator(Elements.find(Col));
      }

      const_iterator find(int Col) const
      {
         return const_iterator(Elements.find(Col));
      }

      void clear() { Elements.clear(); }

      template<typename... Args>
      void emplace(int Col, Args&&... args)
      {
         bool Inserted = Elements.emplace(std::forward<Args>(args)...).second;
         CHECK(Inserted);
      }

      void insert(int Col, T const& value)
      {
         bool Inserted = Elements.insert(std::make_pair(Col, value)).second;
         CHECK(Inserted);
      }

      void insert(int Col, T&& value)
      {
         bool Inserted = Elements.insert(std::make_pair(Col, std::move(value))).second;
         CHECK(Inserted);
      }

      void add(int Col, T const& value)
      {
         base_iterator I = Elements.find(Col);
         if (I == Elements.end())
         {
            Elements.insert(std::make_pair(Col, copy(value)));
         }
         else
         {
            I->second += value;
         }
      }

      template <typename U>
      void add(int Col, U&& value)
      {
         base_iterator I = Elements.find(Col);
         if (I == Elements.end())
         {
            Elements.insert(std::make_pair(Col, T(std::move(value))));
         }
         else
         {
            I->second += std::move(value);
         }
      }

      void subtract(int Col, T const& value)
      {
         base_iterator I = Elements.find(Col);
         if (I == Elements.end())
         {
            Elements.insert(std::make_pair(Col, -value));
         }
         else
         {
            I->second -= value;
         }
      }

      template <typename U>
      void subtract(int Col, U&& value)
      {
         base_iterator I = Elements.find(Col);
         if (I == Elements.end())
         {
            Elements.insert(std::make_pair(Col, -std::move(value)));
         }
         else
         {
            I->second -= std::move(value);
         }
      }

      iterator at(int Col) { return iterator(Elements.find(Col)); }
      const_iterator at(int Col) const { return iterator(Elements.find(Col)); }

   private:
      std::map<int, T> Elements;
};

template <typename T, typename U, typename V>
BasicSparseVector<T, U, V>& operator+=(BasicSparseVector<T, U, V>& x, BasicSparseVector<T, U, V> const& Other)
{
   for (auto const& e : Other)
   {
      x.add(e.index, e.value);
   }
   return x;
}

template <typename T, typename U, typename V>
BasicSparseVector<T, U, V>& operator+=(BasicSparseVector<T, U, V>& x, BasicSparseVector<T, U, V>&& Other)
{
   for (auto&& e : Other)
   {
      x.add(e.index, std::move(e.value));
   }
   Other.clear();
   return x;
}

template <typename T, typename U, typename V, typename W, typename X, typename R>
void
inner_prod(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, R&& result)
{
   bool Set = false;
   auto xi = x.begin();
   auto yi = y.begin();
   while (xi != x.end() && yi != y.end())
   {
      if (xi.index() == yi.index())
      {
	 if (!Set)
	 {
	    inner_prod(xi.value(), yi.value(), static_cast<R&&>(result));
	    Set = true;
	 }
	 else
	    add_inner_prod(xi.value(), yi.value(), static_cast<R&&>(result));
         ++xi;
         ++yi;
      }
      // can't allow fall through here since we need to loop back and check against end()
      else if (xi.index() < yi.index())
         ++xi;
      else // if (yi.index < xi.index)
         ++yi;
   }
   if (!Set)
   {
      // TODO: implement this
      //      clear(result);
      PANIC("Not implemented");
   }
}

template <typename T, typename U, typename V, typename W, typename X, typename R>
void
add_inner_prod(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, R&& result)
{
   auto xi = x.begin();
   auto yi = y.begin();
   while (xi != x.end() && yi != y.end())
   {
      if (xi.index() == yi.index())
      {
         add_inner_prod(xi.value(), yi.value(), static_cast<R&&>(result));
         ++xi;
         ++yi;
      }
      // can't allow fall through here since we need to loop back and check against end()
      else if (xi.index() < yi.index())
         ++xi;
      else // if (yi.index < xi.index)
         ++yi;
   }
}

template <typename T>
using SparseRowVector = BasicSparseVector<T, detail::SparseRowIterator<T>, detail::ConstSparseRowIterator<T>>;

// A row of a sparse matrix.  This attaches the row number to a sparse vector.
// The row number is fixed at construction and is thereafter immutable.  hence
// most operations on the row will delegate to the base class, slicing it if
// it is an l-value, which is harmless in this context.
template <typename T>
class SparseMatrixRow : public SparseRowVector<T>
{
   public:
      SparseMatrixRow() = delete;
      explicit SparseMatrixRow(int Row_) : Row(Row_) {}

      int row() const { return Row; }

   private:
      int Row;
};

template <typename T>
class SparseMatrix
{
   public:
      using value_type     = T;
      using row_type       = SparseMatrixRow<T>;
      using iterator       = typename std::vector<SparseMatrixRow<T>>::iterator;
      using const_iterator = typename std::vector<SparseMatrixRow<T>>::const_iterator;

      SparseMatrix() = delete;

      SparseMatrix(SparseMatrix&& other) = default;

      SparseMatrix(SparseMatrix const&) = delete;

      SparseMatrix(int Rows_, int Cols_) : Cols(Cols_)
      {
         RowStorage.reserve(Rows_);
         for (int i = 0; i < Rows_; ++i)
         {
            RowStorage.emplace_back(i);
         }
      }

      SparseMatrix& operator=(SparseMatrix<T>&& x)
      {
	 Cols = x.Cols;
	 RowStorage = std::move(x.RowStorage);
	 x.Cols = 0;
	 return *this;
      }

      int rows() const { return RowStorage.size(); }
      int cols() const { return Cols; }
      int nnz() const { int Result = 0; for (auto const& x : RowStorage) { Result += x.nnz(); } return Result; }

      void clear() { for (auto& r : RowStorage) { r.clear(); } }

      bool empty() const { return RowStorage.empty(); }

      iterator begin() { return RowStorage.begin(); }
      iterator end() { return RowStorage.end(); }

      const_iterator begin() const { return RowStorage.begin(); }
      const_iterator end() const { return RowStorage.end(); }

      const_iterator cbegin() const { return RowStorage.begin(); }
      const_iterator cend() const { return RowStorage.end(); }

      row_type& operator[](int r) { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }
      row_type const& operator[](int r) const { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }

      row_type& row(int r) { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }
      row_type const& row(int r) const { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }

      template<typename... Args>
      void emplace(int Row, int Col, Args&&... args)
      {
         bool Inserted = RowStorage[Row].emplace(std::forward<Args>(args)...).second;
         CHECK(Inserted);
      }

      template <typename U>
      void insert(int r, int c, U const& value)
      {
         DEBUG_RANGE_CHECK_OPEN(r, 0, RowStorage.size());
         DEBUG_RANGE_CHECK_OPEN(c, 0, Cols);
         RowStorage[r].insert(c, value);
      }

      template <typename U>
      void insert(int r, int c, U&& value)
      {
         DEBUG_RANGE_CHECK_OPEN(r, 0, RowStorage.size());
         DEBUG_RANGE_CHECK_OPEN(c, 0, Cols);
         RowStorage[r].insert(c, std::move(value));
      }

      template <typename U>
      void add(int r, int c, U const& value)
      {
         DEBUG_RANGE_CHECK_OPEN(r, 0, RowStorage.size());
         DEBUG_RANGE_CHECK_OPEN(c, 0, Cols);
         RowStorage[r].add(c, value);
      }

      template <typename U>
      void add(int r, int c, U&& value)
      {
         DEBUG_RANGE_CHECK(r, 0, RowStorage.size());
         DEBUG_RANGE_CHECK(c, 0, Cols);
         RowStorage[r].add(c, std::move(value));
      }

      template <typename U>
      void subtract(int r, int c, U const& value)
      {
         DEBUG_RANGE_CHECK_OPEN(r, 0, RowStorage.size());
         DEBUG_RANGE_CHECK_OPEN(c, 0, Cols);
         RowStorage[r].subtract(c, value);
      }


      template <typename U>
      void subtract(int r, int c, U&& value)
      {
         DEBUG_RANGE_CHECK_OPEN(r, 0, RowStorage.size());
         DEBUG_RANGE_CHECK_OPEN(c, 0, Cols);
         RowStorage[r].subtract(c, std::move(value));
      }

      template <typename U>
      SparseMatrix operator*=(U const& x)
      {
	 for (auto& r : RowStorage)
	 {
	    for (auto&& c : r)
	    {
	       c.value *= x;
	    }
	 }
	 return *this;
      }

      std::vector<SparseMatrixRow<T>> const& access_rows() const { return RowStorage; }
      std::vector<SparseMatrixRow<T>>& access_rows() { return RowStorage; }

   private:
      int Cols;
      std::vector<SparseMatrixRow<T>> RowStorage;
};

template <typename T>
SparseMatrix<T>& operator+=(SparseMatrix<T>& x, SparseMatrix<T> const& Other)
{
   CHECK_EQUAL(x.rows(), Other.rows());
   CHECK_EQUAL(x.cols(), Other.cols());

   int Sz = x.rows();
   for (int r = 0; r < Sz; ++r)
   {
      x[r] += Other[r];
   }
   return x;
}

template <typename T>
SparseMatrix<T>& operator+=(SparseMatrix<T>& x, SparseMatrix<T>&& Other)
{
   CHECK_EQUAL(x.rows(), Other.rows());
   CHECK_EQUAL(x.cols(), Other.cols());

   int Sz = x.rows();
   for (int r = 0; r < Sz; ++r)
   {
      x[r] += std::move(Other[r]);
   }
   Other.clear();
   return x;
}

template <typename T>
SparseMatrix<T>& operator*=(SparseMatrix<T>& x, T const& a)
{
   for (auto& r : x)
   {
      for (auto& c : r)
      {
         c *= a;
      }
   }
   return x;
}

template <typename T>
inline
std::enable_if_t<std::is_copy_constructible<T>::value, T>
copy(T const& x)
{
   return x;
}

// for a temporary, the version that takes a universal reference is found first
// in overload resolution and we can avoid the copy im this case
template <typename T>
inline
std::enable_if_t<std::is_move_constructible<T>::value && !std::is_lvalue_reference<T>::value, T>
copy(T&& x) // x is a universal reference here, which isn't really what we want
{
   return std::move(x);
}

template <typename T>
inline
SparseMatrix<T>
copy(SparseMatrix<T> const& x)
{
   SparseMatrix<T> Result(x.rows(), x.cols());
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
	 Result.insert(r.row(), c.col(), copy(c.value));
      }
   }
   return Result;
}

// A proxy class to represent the Hermitian conjugate of a matrix
template <typename T>
struct SparseMatrixHermitian
{
   using base_type  = T;
   using value_type = typename T::value_type;
   using reference  = T const&;

   explicit SparseMatrixHermitian(T const& x) : x_(x) {}

   reference base() const { return x_; }

   private:
      reference x_;
};

template <typename T>
inline
SparseMatrixHermitian<SparseMatrix<T>>
herm(SparseMatrix<T> const& x)
{
   return SparseMatrixHermitian<SparseMatrix<T>>(x);
}

template <typename T, typename U>
SparseMatrix<remove_proxy_t<decltype(std::declval<T>() * std::declval<U>())>>
operator*(SparseMatrix<T> const& x, SparseMatrix<U> const& y)
{
   SparseMatrix<remove_proxy_t<decltype(std::declval<T>() * std::declval<U>())>> Result(x.rows(), y.cols());
   for (auto const& rx : x)
   {
      for (auto const& cx : rx)
      {
	 for (auto const& cy : y.row(cx.col()))
	 {
	    Result.add(rx.row(), cy.col(), cx.value * cy.value);
	 }
      }
   }
   return Result;
}

template <typename T, typename U>
SparseMatrix<remove_proxy_t<decltype(std::declval<T>() * std::declval<U>())>>
operator*(SparseMatrix<T> const& x, SparseMatrixHermitian<SparseMatrix<U>> const& y)
{
   using element_type = remove_proxy_t<decltype(std::declval<T>() * std::declval<U>())>;
   SparseMatrix<element_type> Result(x.rows(), y.base().rows());
   for (auto const& rx : x)
   {
      for (auto const& ry : y.base())
      {
	 bool HasElement = false;
	 auto xi = rx.begin();
	 auto yi = ry.begin();
	 // find the first occurance of an element with xi.col == yi.col
	 while (xi != rx.end() && yi != ry.end() && (xi.col() != yi.col()))
	 {
	    if (xi.col() < yi.col())
	       ++xi;
	    else if (yi.col() < xi.col())
	       ++yi;
	 }
	 if (xi != rx.end() && yi != ry.end())
	 {
	    // we have an element
	    element_type E = rx.value * ry.value;
	    // search for others
	    while (xi != rx.end() && yi != ry.end())
	    {
	       if (xi.col() == yi.col())
	       {
		  E += rx.value * ry.value;
		  ++xi;
		  ++yi;
	       }
	       else if (xi.col() < yi.col())
		  ++xi;
	       else if (yi.col() < xi.col())
		  ++yi;
	    }
	    Result.insert(rx.row(), ry.row(), std::move(E));
	 }
      }
   }
   return Result;
}

template <typename T, typename U>
SparseMatrix<remove_proxy_t<decltype(std::declval<T>() * std::declval<U>())>>
operator*(SparseMatrixHermitian<SparseMatrix<U>> const& x, SparseMatrix<T> const& y)
{
   using element_type = remove_proxy_t<decltype(std::declval<T>() * std::declval<U>())>;
   SparseMatrix<element_type> Result(x.rows(), y.base().rows());
   for (auto const& rx : x.base())
   {
      for (auto const& cx : rx)
      {
	 for (auto const& cy : y.row(cx.col()))
	 {
	    Result.add(rx.row(), cy.col(), cx.value * cy.value);
	 }
      }
   }
   return Result;
}

inline
void write_matrix_market(std::ostream& out, SparseMatrix<double> const& x)
{
   out << x.rows() << ' ' << x.cols() << ' ' << x.nnz() << '\n';
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
         out << (r.row()+1) << ' ' << (c.col()+1) << ' ' << c << '\n';
      }
   }
}

template <typename T>
inline
std::ostream& operator<<(std::ostream& out, SparseMatrix<T> const& x)
{
   out << "sparse matrix [" << x.rows() << ' ' << x.cols() << ' ' << x.nnz() << "]\n";
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
         out << "element at " << r.row() << ' ' << c.col() << ":\n"
             << c.value << '\n';
      }
   }
   return out;
}

inline
std::ostream& operator<<(std::ostream& out, SparseMatrix<double> const& x)
{
   out << '[' << x.rows() << ' ' << x.cols() << ' ' << x.nnz() << "]\n";
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
         out << r.row() << ' ' << c.col() << ' ' << c << '\n';
      }
   }
   return out;
}

} // namespace blas

template <typename T>
struct ScalarTypes<blas::SparseMatrix<T>> : ScalarTypes<T> {};

#endif
