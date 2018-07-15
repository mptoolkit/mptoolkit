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
#include "functors.h"
#include "diagonalmatrix.h"
#include "number_traits.h"
#include <vector>
#include <map>
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#include "pstreamio.h"
#endif

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

      static_assert(std::is_nothrow_move_constructible<value_type>::value, "");

      BasicSparseVector() = default;
      BasicSparseVector(BasicSparseVector&& Other) noexcept = default;
      BasicSparseVector& operator=(BasicSparseVector&& Other) = default;
      BasicSparseVector& operator=(BasicSparseVector const& Other) = default;
      ~BasicSparseVector() = default;

      BasicSparseVector(BasicSparseVector const&) = default;

#if 0
      // this is bogus?
      template <typename U>
      explicit BasicSparseVector(BasicSparseVector const& Other)
	 : Elements(Other.Elements) {}
#endif

      int nnz() const { return Elements.size(); }

      bool empty() const { return Elements.empty(); }
      bool is_zero() const { return Elements.empty(); }

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
         bool Inserted = Elements.emplace(std::piecewise_construct, std::forward_as_tuple(Col),
					  std::forward_as_tuple(std::forward<Args>(args)...)).second;
         DEBUG_CHECK(Inserted);
      }

      template <typename U>
      iterator insert(int Col, U&& value)
      {
	 auto i = Elements.insert(std::make_pair(Col, std::forward<U>(value)));
         //bool Inserted = Elements.insert(std::make_pair(Col, std::forward<U>(value))).second;
         DEBUG_CHECK(i.second);
	 return iterator(i.first);
      }

      void add(int Col, T const& value)
      {
         base_iterator I = Elements.find(Col);
         if (I == Elements.end())
         {
	    using ::copy;
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
            Elements.emplace(std::make_pair(Col, T(std::move(value))));
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

      bool exists(int Col) const { return Elements.find(Col) != Elements.end(); }

      void erase(int Col) { Elements.erase(Col); }

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

template <typename T, typename U, typename V>
BasicSparseVector<T, U, V>& operator-=(BasicSparseVector<T, U, V>& x, BasicSparseVector<T, U, V> const& Other)
{
   for (auto const& e : Other)
   {
      x.add(e.index, -e.value);
   }
   return x;
}

template <typename T, typename U, typename V>
BasicSparseVector<T, U, V>& operator-=(BasicSparseVector<T, U, V>& x, BasicSparseVector<T, U, V>&& Other)
{
   for (auto&& e : Other)
   {
      x.add(e.index, -std::move(e.value));
   }
   Other.clear();
   return x;
}

template <typename T, typename U, typename V, typename W, typename X, typename R, typename Nested, typename Tag>
auto
inner_prod_nested_tag(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, R&& Result, Nested&& Nest, Tag)
{
   // device version, allocate a device vector to accumulate the partial results
   using result_value = remove_proxy_t<R>;

   // determine the maximum size we need
   auto xlast = x.end();
   auto ylast = y.end();
   if (xlast == x.begin() || ylast == y.begin())
   {
      clear(Result);
      return;
   }
   xlast--;
   ylast--;
   auto Size = std::max(xlast.index(), ylast.index());
   blas::Vector<result_value, Tag> Acc(Size, 0.0);
   bool Set = false;
   auto xi = x.begin();
   auto yi = y.begin();
   while (xi != x.end() && yi != y.end())
   {
      if (xi.index() == yi.index())
      {
	 Nest(xi.value(), yi.value(), Acc[xi.index()]);
         ++xi;
         ++yi;
      }
      // can't allow fall through here since we need to loop back and check against end()
      else if (xi.index() < yi.index())
         ++xi;
      else // if (yi.index < xi.index)
         ++yi;
   }
   sum(Acc, Result);
}

template <typename T, typename U, typename V, typename W, typename X, typename R, typename Nested>
void
inner_prod_nested_tag(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, R&& Result,
		      Nested&& Nest, cpu_tag)
{
   // CPU version, accumulate results into a value
   using result_value = remove_proxy_t<decltype(std::declval<Nested>()(std::declval<T>(),std::declval<T>()))>;
   bool Set = false;
   auto xi = x.begin();
   auto yi = y.begin();
   while (xi != x.end() && yi != y.end())
   {
      if (xi.index() == yi.index())
      {
	 if (!Set)
	 {
	    Result = Nest(xi.value(), yi.value());
	    Set = true;
	 }
	 else
	    Result += Nest(xi.value(), yi.value());
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

template <typename T, typename U, typename V, typename W, typename X, typename R, typename Nested>
auto
inner_prod_nested(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, R&& Result, Nested&& Nest)
{
   return inner_prod_nested_tag(x, y, std::forward<R>(Result), std::forward<Nested>(Nest), tag_of_t<T>());
}

// The case where T is a device type but we want the result on the CPU
template <typename T, typename U, typename V, typename W, typename X, typename Nested, typename Tag>
auto
inner_prod_nested_tag(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, Nested&& Nest, Tag)
{
   // device version, allocate a device vector to accumulate the partial results
   using result_value = remove_proxy_t<decltype(std::declval<Nested>()(std::declval<T>(),std::declval<T>()))>;

   // determine the maximum size we need
   auto xlast = x.end();
   auto ylast = y.end();
   if (xlast == x.begin() || ylast == y.begin())
   {
      return result_value{};
   }
   xlast--;
   ylast--;
   auto Size = std::max(xlast.index(), ylast.index());
   blas::Vector<result_value, Tag> Acc(Size, 0.0);
   bool Set = false;
   auto xi = x.begin();
   auto yi = y.begin();
   while (xi != x.end() && yi != y.end())
   {
      if (xi.index() == yi.index())
      {
	 Nest(xi.value(), yi.value(), Acc[xi.index()]);
         ++xi;
         ++yi;
      }
      // can't allow fall through here since we need to loop back and check against end()
      else if (xi.index() < yi.index())
         ++xi;
      else // if (yi.index < xi.index)
         ++yi;
   }
   return get_wait(sum(Acc));
}

template <typename T, typename U, typename V, typename W, typename X, typename Nested>
inline
auto
inner_prod_nested_tag(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, Nested&& Nest, cpu_tag)
{
   using result_type = remove_proxy_t<decltype(std::declval<Nested>()(std::declval<T>(),std::declval<T>()))>;
   result_type Acc(0.0);
   inner_prod_nested(x, y, Acc, std::forward<Nested>(Nest));
   return Acc;
}

template <typename T, typename U, typename V, typename W, typename X, typename Nested>
inline
auto
inner_prod_nested(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, Nested&& Nest)
{
   return inner_prod_nested_tag(x, y, std::forward<Nested>(Nest), tag_of_t<T>());
}

template <typename T, typename U, typename V, typename W, typename X, typename R, typename Nested>
void
add_inner_prod_nested(BasicSparseVector<T,U,V> const& x, BasicSparseVector<T,W,X> const& y, R&& result, Nested&& Nest)
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
      using row_reference  = row_type&;
      using const_row_reference = row_type const&;
      using iterator       = typename std::vector<SparseMatrixRow<T>>::iterator;
      using const_iterator = typename std::vector<SparseMatrixRow<T>>::const_iterator;

      SparseMatrix() noexcept  : Cols(0) {}

      SparseMatrix(SparseMatrix&& other) noexcept = default;

      SparseMatrix(SparseMatrix const&) = delete;

      template <typename U>
      explicit SparseMatrix(SparseMatrix<U> const& Other)
	 : Cols(Other.Cols), RowStorage(Other.RowStorage) {}

      template <typename U>
      explicit SparseMatrix(DiagonalMatrix<U> const& Other);

      explicit SparseMatrix(DiagonalMatrix<T>&& Other);

      SparseMatrix(int Rows_, int Cols_) : Cols(Cols_)
      {
         RowStorage.reserve(Rows_);
         for (int i = 0; i < Rows_; ++i)
         {
            RowStorage.emplace_back(i);
         }
      }

      SparseMatrix& operator=(SparseMatrix<T>&& x) noexcept
      {
	 Cols = x.Cols;
	 RowStorage = std::move(x.RowStorage);
	 x.Cols = 0;
	 return *this;
      }

      SparseMatrix& operator=(SparseMatrix<T> const& x) = default;

      int rows() const { return RowStorage.size(); }
      int cols() const { return Cols; }
      int nnz() const { int Result = 0; for (auto const& x : RowStorage) { Result += x.nnz(); } return Result; }

      void clear() { for (auto& r : RowStorage) { r.clear(); } }

      bool empty() const { return RowStorage.empty(); }

      bool is_zero() const
      {
	 for (auto const& r : RowStorage)
	 {
	    if (!r.empty())
	       return false;
	 }
	 return true;
      }

      void erase(int r, int c)
      {
	 DEBUG_RANGE_CHECK_OPEN(r,0,RowStorage.size());
	 DEBUG_RANGE_CHECK_OPEN(c,0,Cols);
	 RowStorage[r].erase(c);
      }

      bool exists(int r, int c) const
      {
	 DEBUG_RANGE_CHECK_OPEN(r,0,RowStorage.size());
	 DEBUG_RANGE_CHECK_OPEN(c,0,Cols);
	 return RowStorage[r].exists(c);
      }

      iterator begin() { return RowStorage.begin(); }
      iterator end() { return RowStorage.end(); }

      const_iterator begin() const { return RowStorage.begin(); }
      const_iterator end() const { return RowStorage.end(); }

      const_iterator cbegin() const { return RowStorage.begin(); }
      const_iterator cend() const { return RowStorage.end(); }

      row_type& operator[](int r) { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }
      row_type const& operator[](int r) const { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }

      row_reference row(int r) { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }
      const_row_reference row(int r) const { DEBUG_RANGE_CHECK(r, 0, RowStorage.size()); return RowStorage[r]; }

      template<typename... Args>
      void emplace(int Row, int Col, Args&&... args)
      {
         RowStorage[Row].emplace(Col, std::forward<Args>(args)...);
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
      SparseMatrix& operator*=(U const& x)
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

// numeric_type_of specialization
template <typename T>
struct numeric_type_of<SparseMatrix<T>> : numeric_type_of<T> {};

#if defined(USE_PSTREAM)
template <int Format, typename T>
PStream::opstreambuf<Format>&
operator<<(PStream::opstreambuf<Format>& out, SparseMatrix<T> const& M)
{
   using st = typename PStream::opstreambuf<Format>::size_type;
   st Rows = M.rows(), Cols = M.cols(), Nnz = M.nnz();
   out << MatrixFormats::Coordinate;
   out << Rows << Cols << Nnz;

   for (auto const& r : M)
   {
      for (auto const& c : r)
      {
	 st row = r.row();
	 st col = c.col();
	 out << row << col;
	 out << c.value;
      }
   }
   return out;
}

template <int Format, typename T>
PStream::ipstreambuf<Format>&
operator>>(PStream::ipstreambuf<Format>& in, SparseMatrix<T>& M)
{
   MatrixFormats::type t;
   in >> t;
   if (t == MatrixFormats::Coordinate)
   {
      using st = typename PStream::opstreambuf<Format>::size_type;
      st Rows, Cols, Nnz;
      in >> Rows >> Cols;
      SparseMatrix<T> Result(Rows, Cols);
      in >> Nnz;
      for (unsigned i = 0; i < Nnz; ++i)
      {
	 st r,c;
	 in >> r >> c;
	 Result.insert(r, c, in.template read<T>());
      }
      M = std::move(Result);
   }
   else
   {
      PANIC("Format not implemented")(t);
   }
   return in;
}

#if 0
operator<<(PStream::opstreambuf<Format>& out, SparseMatrix<T> const& M)
{
   typedef typename PStream::opstreambuf<Format>::size_type st;
   st Rows = M.rows(), Cols = M.cols(), Nnz = M.nnz();
   out << MatrixFormats::Coordinate;
   out << Rows << Cols << Nnz;
   for (auto const& I : M)
   {
      for (auto const& J : I)
      {
	 st i = I.row();
	 st j = J.col();
	 out << i << j << J.value;
      }
   }
   return out;
}

template <typename T>
struct stream_extract<blas::SparseMatrix<T>>
{
   template <typename Format>
   blas::SparseMatrix<T> operator()(PStream::ipstreambuf<Format>& in) const
   {
      typedef typename PStream::ipstreambuf<Format>::size_type st;
      st Rows, Cols, Nnz;
      blas::MatrixFormats::type t;
      in >> t;
      if (t != blas::MatrixFormats::Coordinate)
      {
	 PANIC("Unsupported sparse matrix format");
      }
      in >> Rows >> Cols >> Nnz;
      blas::SparseMatrix<T> Result(Rows, Cols);
      for (unsigned n = 0; n < Nnz; ++n)
      {
	 st i, j;
	 in >> i >> j;
	 Result.insert(i,j, PStream::extract<T>(in));
      }
      return Result;
   }
};
#endif

#endif

template <typename T, typename U>
SparseMatrix<T>& operator+=(SparseMatrix<T>& x, SparseMatrix<U> const& Other)
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

template <typename T, typename U>
SparseMatrix<T>& operator+=(SparseMatrix<T>& x, SparseMatrix<U>&& Other)
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

template <typename T, typename U>
SparseMatrix<T>& operator-=(SparseMatrix<T>& x, SparseMatrix<U> const& Other)
{
   CHECK_EQUAL(x.rows(), Other.rows());
   CHECK_EQUAL(x.cols(), Other.cols());

   int Sz = x.rows();
   for (int r = 0; r < Sz; ++r)
   {
      x[r] -= Other[r];
   }
   return x;
}

template <typename T, typename U>
SparseMatrix<T>& operator-=(SparseMatrix<T>& x, SparseMatrix<U>&& Other)
{
   CHECK_EQUAL(x.rows(), Other.rows());
   CHECK_EQUAL(x.cols(), Other.cols());

   int Sz = x.rows();
   for (int r = 0; r < Sz; ++r)
   {
      x[r] -= std::move(Other[r]);
   }
   Other.clear();
   return x;
}


template <typename T>
SparseMatrix<T> operator-(SparseMatrix<T> const& x)
{
   SparseMatrix<T> Result(x.rows(), x.cols());
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
	 Result.insert(r.row(), c.col(), -c.value);
      }
   }
   return Result;
}

template <typename T>
SparseMatrix<T> operator-(SparseMatrix<T>&& x)
{
   for (auto& r : x)
   {
      for (auto&& c : r)
      {
	 c.value = -std::move(c.value);
      }
   }
   return std::move(x);
}

template <typename T, typename U>
SparseMatrix<remove_proxy_t<decltype(std::declval<T>() + std::declval<U>())>>
operator+(SparseMatrix<T> const& x, SparseMatrix<U> const& y)
{
   using ::copy;
   SparseMatrix<remove_proxy_t<decltype(std::declval<T>() + std::declval<U>())>>
      Result(copy(x));
   Result += y;
   return Result;
}

template <typename T, typename U>
SparseMatrix<remove_proxy_t<decltype(std::declval<T>() - std::declval<U>())>>
operator-(SparseMatrix<T> const& x, SparseMatrix<U> const& y)
{
   using ::copy;
   SparseMatrix<remove_proxy_t<decltype(std::declval<T>() + std::declval<U>())>>
      Result(copy(x));
   Result -= y;
   return Result;
}

template <typename T>
inline
SparseMatrix<T>
copy(SparseMatrix<T> const& x)
{
   using ::copy;
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

// conj

template <typename T>
inline
void
inplace_conj(SparseMatrix<T>& x)
{
   for (auto& r : x)
   {
      for (auto&& c : r)
      {
	 inplace_conj(c.value);
      }
   }
}

template <typename T>
inline
SparseMatrix<T>
conj(SparseMatrix<T> const& x)
{
   using std::conj;
   SparseMatrix<T> Result(x.rows(), x.cols());
   for (auto const& r : x)
   {
      for (auto const& c : r)
      {
	 Result.insert(r.row(), c.col(), conj(c.value));
      }
   }
   return Result;
}

template <typename T>
inline
SparseMatrix<T>
conj(SparseMatrix<T>&& x)
{
   inplace_conj(x);
   return std::move(x);
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
operator*(SparseMatrix<T> const& x, DiagonalMatrix<U> const& y)
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
operator*(DiagonalMatrix<T> const& x, SparseMatrix<U> const& y)
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
	    element_type E = xi.value() * yi.value();
	    // search for others
	    while (xi != rx.end() && yi != ry.end())
	    {
	       if (xi.col() == yi.col())
	       {
		  E += xi.value() * yi.value();
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
   SparseMatrix<element_type> Result(x.base().cols(), y.cols());
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

template <typename T>
template <typename U>
SparseMatrix<T>::SparseMatrix(DiagonalMatrix<U> const& Other)
   : Cols(Other.cols()), RowStorage(Other.rows())
{
   int c = 0;
   for (auto const& x : Other.diagonal())
   {
      this->insert(c, c, copy(x));
      ++c;
   }
}

template <typename T>
SparseMatrix<T>::SparseMatrix(DiagonalMatrix<T>&& Other)
   : Cols(Other.cols()), RowStorage(Other.rows())
{
   int c = 0;
   for (auto&& x : Other.diagonal())
   {
      this->insert(c, c, copy(x));
      ++c;
   }
}

} // namespace blas

template <typename T>
struct ScalarTypes<blas::SparseMatrix<T>> : ScalarTypes<T> {};

#endif
