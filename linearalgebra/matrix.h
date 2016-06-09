// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrix.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MATRIX_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define MATRIX_H_DHJFKHURWEIRY4572893475489YRUI34897

#include "matrixfwd.h"
#include "matrixinterface.h"
#include "matrixoperations.h"
#include "datablockreference.h"
#include "matrixptriterator.h"
#include "crtp_matrix.h"
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

namespace LinearAlgebra
{

struct MatrixDimensions
{
   size_type size1, size2;
   difference_type stride1, stride2;

   MatrixDimensions() {}

   // There is a corner case for Nx1 matrices (column vectors),
   // where we could have stride1 == stride2 == 1, but in this case
   // transpose is trivial, and the transpose must have a stride1 that
   // is not 1.  Thus, we do a check and set the original stride2 to be size1 for this case.
   // This should be otherwise harmless.
   MatrixDimensions(size_type Rows, size_type Cols, RowMajor) 
      : size1(Rows), size2(Cols), stride1(Cols), stride2(Cols == 1 ? Rows : 1) {}
   //      : size1(Rows), size2(Cols), stride1(Cols), stride2(Rows) {}

   // FIXME: is this correct?
   MatrixDimensions(size_type Rows, size_type Cols, ColMajor) 
      : size1(Rows), size2(Cols), stride1(Rows == 1 ? Cols : 1), stride2(Rows) {}
   //      : size1(Rows), size2(Cols), stride1(Cols), stride2(Rows) {}
   
   // TODO: the above note applies to this ctor too.
   //   MatrixDimensions(size_type Rows, size_type Cols, 
   //		    difference_type RowStride, difference_type ColStride)
   //      : size1(Rows), size2(Cols), stride1(RowStride), stride2(ColStride) {}

   size_type size() const { return size1*size2; }
};

#if defined(USE_PSTREAM)
template <int Format>
PStream::opstreambuf<Format>& 
operator<<(PStream::opstreambuf<Format>& out, MatrixDimensions const& d);

template <int Format>
PStream::ipstreambuf<Format>& 
operator>>(PStream::opstreambuf<Format>& in, MatrixDimensions& d);
#endif

template <typename Scalar, typename Orientation = RowMajor, typename Derived = void>
class MatrixRef;

template <typename Scalar, typename Orient, typename Derived>
struct MatrixRefDerivedType
{
   typedef Derived type;
};

template <typename Scalar, typename Orient>
struct MatrixRefDerivedType<Scalar, Orient, void>
{
   typedef MatrixRef<Scalar, Orient> type;
};

template <typename Scalar, typename Orientation, typename Derived>
class MatrixRef : public MatrixBase<typename MatrixRefDerivedType<Scalar, Orientation, Derived>::type>
//MatrixBase<MatrixRef<Scalar, Orientation> >
{
   private:
      typedef BlockReference<Scalar, NoHeader, MatrixDimensions> block_type;

      typedef MatrixBase<typename MatrixRefDerivedType<Scalar, Orientation, Derived>::type> base_type;

   public:
      using base_type::operator();

      typedef Scalar value_type;
      typedef Scalar& reference;
      typedef Scalar* pointer;

      typedef Scalar const& const_reference;
      typedef Scalar const* const_pointer;

      MatrixRef(MatrixRef const& V) 
	: Block(V.Block) {}

      MatrixRef& operator=(MatrixRef const& V)
      {
	 this->assign(V);   // safe because assign() calls begin(), which calls cow()
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixRef&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_matrix<U>, MatrixRef&>::type
      operator=(NoAliasProxy<U> const& x);

      size_type size() const { return Block.size(); }
      size_type size1() const { return Block.local_header().size1; }
      size_type size2() const { return Block.local_header().size2; }

      Scalar const& operator()(size_type i, size_type j) const 
      { DEBUG_RANGE_CHECK_OPEN(i, 0U, this->size1()); 
        DEBUG_RANGE_CHECK_OPEN(j, 0U, this->size2());
	return *(Block.get() + this->stride1() * i + this->stride2() * j); }

      Scalar& operator()(size_type i, size_type j)
      { DEBUG_RANGE_CHECK_OPEN(i, 0U, this->size1()); 
        DEBUG_RANGE_CHECK_OPEN(j, 0U, this->size2());
        this->cow(); return *(Block.get() + this->stride1() * i + this->stride2() * j); }

      // members defined in MatrixStride
      Scalar* data() { this->cow(); return Block.get(); }
      Scalar const* data() const { return Block.get(); }
      difference_type stride1() const { return Block.local_header().stride1; }
      difference_type stride2() const { return Block.local_header().stride2; }

      // members defined in MatrixRef
      void inplace_transpose();

      difference_type leading_dimension() const 
         { return std::max(Block.local_header().stride1, Block.local_header().stride2); }

   protected:
      block_type const& block() const { return Block; } 
      block_type& block() { return Block; }  // danger: this does not call cow()

      MatrixRef() : Block(NoHeader(), MatrixDimensions(0,0, Orientation())) {}

      explicit MatrixRef(size_type Rows, size_type Cols)
	: Block(Rows*Cols, NoHeader(), MatrixDimensions(Rows, Cols, Orientation())) {}

      explicit MatrixRef(size_type Rows, size_type Cols, Scalar const& Fill)
                  : Block(Rows*Cols, Fill, NoHeader(), MatrixDimensions(Rows, Cols, Orientation())) 
               { }
         //         : Block(Rows*Cols, NoHeader(), MatrixDimensions(Rows, Cols, Orientation())) 
   //{ LinearAlgebra::Fill<MatrixRef<Scalar, Orientation, Derived>&, Scalar>()(*this, Fill); }
  
      explicit MatrixRef(block_type const& Block_) : Block(Block_) {}

      void resize(size_type NewRows, size_type NewCols);
  
      void swap(MatrixRef& Other) 
        { this->cow(); Other.cow(); std::swap(Block, Other.Block); }

   private:
      void cow() { Block.cow(); }

      block_type Block;

   template <typename U, typename Orient, typename D> friend class MatrixRef;
   //   template <typename U, typename Orient, > friend class MatrixConstRef;
};

// interface

template <typename Scalar, typename Orientation>
struct interface<MatrixRef<Scalar, Orientation> >
{
   typedef Concepts::ContiguousMatrix<Scalar, Orientation, void> type;
   typedef Scalar value_type;
};

// iterators

template <typename Scalar, typename Orient, typename Derived>
struct Iterate<MatrixRef<Scalar, Orient, Derived>&>
{
   typedef MatrixPtrIterator<Scalar, Orient> result_type;
   typedef MatrixRef<Scalar, Orient, Derived>& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

template <typename Scalar, typename Orient, typename Derived>
struct Iterate<MatrixRef<Scalar, Orient, Derived> >
{
   typedef MatrixPtrIterator<Scalar const, Orient> result_type;
   typedef MatrixRef<Scalar, Orient, Derived> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

//
// Matrix
//
// matrix with copy-on-write semantics.  This is the 'default' matrix class.
//

template <typename Scalar, typename Orientation>
class Matrix : public MatrixRef<Scalar, Orientation, Matrix<Scalar, Orientation> >
{
   public:
      typedef Scalar*       iterator;
      typedef Scalar const* const_iterator;

      Matrix() {}

      Matrix(size_type Rows, size_type Cols);

      Matrix(size_type Rows, size_type Cols, Scalar const& fill_);

      template <typename Iter>
      Matrix(size_type Rows, size_type Cols, Iter first, Iter last);

      Matrix(Matrix const& V)
	: MatrixRef<Scalar, Orientation, Matrix<Scalar, Orientation> >(V.block().copy()) { }

      template <typename U>
      Matrix(U const& x, typename boost::enable_if<is_matrix<U> >::type* dummy = 0);

      Matrix& operator=(Matrix const& V);

      template <typename U>
      typename boost::enable_if<is_matrix<U>, Matrix&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_matrix<U>, Matrix&>::type
      operator=(NoAliasProxy<U> const& x);

      void resize(size_type NewRows, size_type NewCols)
	{ this->MatrixRef<Scalar, Orientation, Matrix<Scalar, Orientation> >::resize(NewRows, NewCols); }

      void swap(Matrix& Other) { this->MatrixRef<Scalar, Orientation>::swap(Other); }
};

// interface

template <typename Scalar, typename Orientation>
struct interface<Matrix<Scalar, Orientation> >
{
   using type = Concepts::ContiguousMatrix<Scalar, Orientation, void>;
   using value_type = Scalar;
};

// iterators

template <typename Scalar, typename Orient>
struct Iterate<Matrix<Scalar, Orient>&>
{
   typedef MatrixPtrIterator<Scalar, Orient> result_type;
   typedef Matrix<Scalar, Orient>& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

template <typename Scalar, typename Orient>
struct Iterate<Matrix<Scalar, Orient> >
{
   typedef MatrixPtrIterator<Scalar const, Orient> result_type;
   typedef Matrix<Scalar, Orient> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return result_type(x.data(), x.size1(), x.stride1(), x.size2(), x.stride2());
   }
};

// resize

template <typename T, typename Orient>
struct Resize<Matrix<T, Orient>&>
{
   typedef void result_type;
   typedef Matrix<T, Orient>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   void operator()(Matrix<T, Orient>& v, size_type r, size_type c) const
   {
      v.resize(r,c);
   }
};

} // namespace LinearAlgebra

#include "matrix.cc"

#if !defined(LINEARALGEBRA_NO_BLAS)
#include "matrixproductblas_impl.cc"
#endif

#endif
