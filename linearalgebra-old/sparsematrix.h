// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/sparsematrix.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(SPARSEMATRIX_H_A636DJBNDJ48976BNM6NM64RJYG86V5)
#define SPARSEMATRIX_H_A636DJBNDJ48976BNM6NM64RJYG86V5

#include "matrixfwd.h"
#include "matrixoperations.h"
#include "compressedmatrix.h"
#include "vector.h"
#include "mapvector.h"        // to construct the 'default' sparse type
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

namespace LinearAlgebra
{

template <typename T, typename InnerType, typename OuterType>
class SparseMatrix<T, RowMajor, InnerType, OuterType>
   : public CompressedMatrixBase<OuterType, SparseMatrix<T, RowMajor, InnerType, OuterType> >
{
   public:
      typedef CompressedMatrixBase<OuterType,
                                   SparseMatrix<T, RowMajor, InnerType, OuterType> > base_type;

      typedef typename base_type::value_type value_type;
      typedef typename base_type::outer_value_type outer_value_type;

      BOOST_MPL_ASSERT((boost::is_same<T, value_type>));

      typedef typename make_reference<OuterType>::type       vec_reference;
      typedef typename make_const_reference<OuterType>::type const_vec_reference;

      SparseMatrix() : InnerSize_(0), Data_() {}

      SparseMatrix(size_type Rows, size_type Cols)
         : InnerSize_(Cols), Data_(Rows, outer_value_type(Cols)) {}

      SparseMatrix(SparseMatrix const& x)
         : InnerSize_(x.inner_size()), Data_(x.Data_) {}

      template <typename U>
      SparseMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type* dummy = 0);

      template <typename U>
      typename Multiply<SparseMatrix<T, RowMajor, InnerType, OuterType>&, U>::result_type
      operator *=(U const& x)
      {
         return Multiply<SparseMatrix<T, RowMajor, InnerType, OuterType>&, U>()(*this, x);
      }

#if 1
      template <typename Arg1, typename Arg2>
      typename MatrixBracket<SparseMatrix, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y) const
      {
         return MatrixBracket<SparseMatrix, Arg1, Arg2>()(*this, x,y);
      }

      template <typename Arg1, typename Arg2>
      typename MatrixBracket<SparseMatrix&, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y)
      {
         return MatrixBracket<SparseMatrix&, Arg1, Arg2>()(*this, x,y);
      }
#endif

      void resize(size_type r, size_type c);

      vec_reference vec() { return Data_; }
      const_vec_reference vec() const { return Data_; }

      size_type inner_size() const { return InnerSize_; }

   private:
      size_type InnerSize_;
      OuterType Data_;
};


// temporary hack
template <typename T>
inline
SparseMatrix<T> identity_matrix(int Dim)
{
   SparseMatrix<T> Ret(Dim, Dim);
   for (int i = 0; i < Dim; ++i)
      Ret(i,i) = 1;
   return Ret;
}

template <typename T, typename InnerType, typename OuterType>
class SparseMatrix<T, ColMajor, InnerType, OuterType>
   : public CompressedMatrixBase<OuterType, SparseMatrix<T, ColMajor, InnerType, OuterType> >
{
   public:
      typedef CompressedMatrixBase<OuterType,
                                   SparseMatrix<T, ColMajor, InnerType, OuterType> > base_type;

      typedef typename base_type::value_type value_type;
      typedef typename base_type::outer_value_type outer_value_type;

      BOOST_MPL_ASSERT((boost::is_same<T, value_type>));

      typedef typename make_reference<OuterType>::type       vec_reference;
      typedef typename make_const_reference<OuterType>::type const_vec_reference;

      SparseMatrix() : InnerSize_(0), Data_() {}

      SparseMatrix(size_type Rows, size_type Cols)
         : InnerSize_(Rows), Data_(Cols, outer_value_type(Rows)) {}

      SparseMatrix(SparseMatrix const& x)
         : InnerSize_(x.inner_size()), Data_(x.Data_) {}

      template <typename U>
      SparseMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type* dummy = 0);

      template <typename U>
      typename Multiply<SparseMatrix<T, ColMajor, InnerType, OuterType>&, U>::result_type
      operator *=(U const& x)
      {
         return Multiply<SparseMatrix<T, ColMajor, InnerType, OuterType>&, U>()(*this, x);
      }

#if 0
      template <typename Arg1, typename Arg2>
      typename MatrixBracket<SparseMatrix, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y) const
      {
         return MatrixBracket<SparseMatrix, Arg1, Arg2>()(*this, x,y);
      }

      template <typename Arg1, typename Arg2>
      typename MatrixBracket<SparseMatrix&, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y)
      {
         return MatrixBracket<SparseMatrix&, Arg1, Arg2>()(*this, x,y);
      }
#endif

      void resize(size_type r, size_type c);

      vec_reference vec() { return Data_; }
      const_vec_reference vec() const { return Data_; }

      size_type inner_size() const { return InnerSize_; }

   private:
      size_type InnerSize_;
      OuterType Data_;
};

template <typename T, typename InnerType, typename OuterType>
template <typename U>
SparseMatrix<T, RowMajor, InnerType, OuterType>::
SparseMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type*)
   : InnerSize_(Size2<U>()(x)), Data_(Size1<U>()(x), outer_value_type(InnerSize_))
{
   assign(*this, x);
}

template <typename T, typename InnerType, typename OuterType>
template <typename U>
SparseMatrix<T, ColMajor, InnerType, OuterType>::
SparseMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type*)
   : InnerSize_(Size1<U>()(x)), Data_(Size2<U>()(x), outer_value_type(InnerSize_))
{
   assign(*this, x);
}

template <typename T, typename InnerType, typename OuterType>
void
SparseMatrix<T, RowMajor, InnerType, OuterType>::resize(size_type r, size_type c)
{
   using LinearAlgebra::resize;
   resize(Data_, r);
   typename iterator<OuterType>::type I = iterate(Data_);
   while (I)
   {
      resize(*I, c);
      ++I;
   }
   InnerSize_ = c;
}

template <typename T, typename InnerType, typename OuterType>
void
SparseMatrix<T, ColMajor, InnerType, OuterType>::resize(size_type r, size_type c)
{
   using LinearAlgebra::resize;
   resize(Data_, c);
   typename iterator<OuterType>::type I = iterate(Data_);
   while (I)
   {
      resize(*I, r);
      ++I;
   }
   InnerSize_ = r;
}

// interface

template <typename T, typename InnerType, typename OuterType, typename Orient>
struct interface<SparseMatrix<T, Orient, InnerType, OuterType> >
{
   typedef T value_type;
   typedef Concepts::CompressedOuterMatrix<value_type, Orient, void> type;
};

// iterators

template <typename T, typename InnerType, typename OuterType, typename Orient>
struct Iterate<SparseMatrix<T, Orient, InnerType, OuterType>&>
{
   typedef typename Iterate<OuterType&>::result_type iter;
   typedef MatrixOuterIterator<iter, Orient> result_type;
   typedef SparseMatrix<T, Orient, InnerType, OuterType>& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.vec()));
   }
};

template <typename T, typename InnerType, typename OuterType, typename Orient>
struct Iterate<SparseMatrix<T, Orient, InnerType, OuterType> >
{
   typedef typename Iterate<OuterType>::result_type iter;
   typedef MatrixOuterIterator<iter, Orient> result_type;
   typedef SparseMatrix<T, Orient, InnerType, OuterType> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.vec()));
   }
};

// iterate_at

template <typename T, typename InnerType, typename OuterType>
struct MatrixIterateAt<SparseMatrix<T, RowMajor, InnerType, OuterType> >
{
   typedef typename const_inner_iterator<SparseMatrix<T, RowMajor, InnerType, OuterType> >::type
      result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename const_iterator<SparseMatrix<T, RowMajor, InnerType, OuterType> >::type
         I = iterate(m);
      I += i;
      return result_type(iterate_at(*I, j), i);
   }
};

template <typename T, typename InnerType, typename OuterType>
struct MatrixIterateAt<SparseMatrix<T, RowMajor, InnerType, OuterType>&>
{
   typedef typename inner_iterator<SparseMatrix<T, RowMajor, InnerType, OuterType> >::type
      result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename iterator<SparseMatrix<T, RowMajor, InnerType, OuterType> >::type
         I = iterate(m);
      I += i;
      return result_type(iterate_at(*I, j), i);
   }
};

template <typename T, typename InnerType, typename OuterType>
struct MatrixIterateAt<SparseMatrix<T, ColMajor, InnerType, OuterType> >
{
   typedef typename const_inner_iterator<SparseMatrix<T, ColMajor, InnerType, OuterType> >::type
      result_type;
   typedef SparseMatrix<T, ColMajor, InnerType, OuterType> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename const_iterator<SparseMatrix<T, ColMajor, InnerType, OuterType> >::type
         I = iterate(m);
      I += j;
      return result_type(iterate_at(*I, i), j);
   }
};

template <typename T, typename InnerType, typename OuterType>
struct MatrixIterateAt<SparseMatrix<T, ColMajor, InnerType, OuterType>&>
{
   typedef typename inner_iterator<SparseMatrix<T, ColMajor, InnerType, OuterType> >::type
      result_type;
   typedef SparseMatrix<T, ColMajor, InnerType, OuterType>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m,
                          second_argument_type i,
                          third_argument_type j) const
   {
      typename iterator<SparseMatrix<T, ColMajor, InnerType, OuterType> >::type
         I = iterate(m);
      I += j;
      return result_type(iterate_at(*I, i), j);
   }
};

// size

template <typename T, typename InnerType, typename OuterType>
struct Size1<SparseMatrix<T, RowMajor, InnerType, OuterType> >
{
   typedef size_type result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType> const& argument_type;
   result_type operator()(argument_type x) const { return size(x.vec()); }
};

template <typename T, typename InnerType, typename OuterType>
struct Size2<SparseMatrix<T, RowMajor, InnerType, OuterType> >
{
   typedef size_type result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return x.inner_size();
   }
};

template <typename T, typename InnerType, typename OuterType>
struct Size1<SparseMatrix<T, ColMajor, InnerType, OuterType> >
{
   typedef size_type result_type;
   typedef SparseMatrix<T, ColMajor, InnerType, OuterType> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return x.inner_size();
   }
};

template <typename T, typename InnerType, typename OuterType>
struct Size2<SparseMatrix<T, ColMajor, InnerType, OuterType> >
{
   typedef size_type result_type;
   typedef SparseMatrix<T, ColMajor, InnerType, OuterType> const& argument_type;
   result_type operator()(argument_type x) const { return size(x.vec()); }
};

// nnz

template <typename T, typename Orient, typename InnerType, typename OuterType>
struct NNZ<SparseMatrix<T, Orient, InnerType, OuterType> >
{
   typedef SparseMatrix<T, Orient, InnerType, OuterType> const& argument_type;
   typedef size_type result_type;

   result_type operator()(argument_type m) const
   {
      size_type n = 0;
      typename const_iterator<SparseMatrix<T, Orient, InnerType, OuterType> >::type
         I = iterate(m);
      while (I)
      {
         n += nnz(*I);
         ++I;
      }
      return n;
   }
};

// is_zero

template <typename T, typename Orient, typename InnerType, typename OuterType>
struct IsZero<SparseMatrix<T, Orient, InnerType, OuterType> >
{
   typedef SparseMatrix<T, Orient, InnerType, OuterType> const& argument_type;
   typedef bool result_type;

   result_type operator()(argument_type m) const
   {
      typename const_iterator<SparseMatrix<T, Orient, InnerType, OuterType> >::type
         I = iterate(m);
      while (I)
      {
         if (!is_zero(*I))
            return false;
         ++I;
      }
      return true;
   }
};

// resize

template <typename T, typename Orient, typename InnerType, typename OuterType>
struct Resize<SparseMatrix<T, Orient, InnerType, OuterType>&>
{
   typedef void result_type;
   typedef SparseMatrix<T, Orient, InnerType, OuterType>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   void operator()(first_argument_type M, size_type r, size_type c) const
   {
      M.resize(r, c);
   }
};

// vector_view

template <typename T, typename InnerType, typename OuterType, typename Orient>
struct VectorView<SparseMatrix<T, Orient, InnerType, OuterType> >
{
   typedef typename make_const_reference<OuterType>::type result_type;
   typedef SparseMatrix<T, Orient, InnerType, OuterType> const& argument_type;
   result_type operator()(argument_type x) const { return x.vec(); }
};

template <typename T, typename InnerType, typename OuterType, typename Orient>
struct VectorView<SparseMatrix<T, Orient, InnerType, OuterType>&>
{
   typedef typename make_reference<OuterType>::type result_type;
   typedef SparseMatrix<T, Orient, InnerType, OuterType>& argument_type;
   result_type operator()(argument_type x) const { return x.vec(); }
};

// set_element

template <typename T, typename InnerType, typename OuterType, typename Value>
struct SetMatrixElement<SparseMatrix<T, RowMajor, InnerType, OuterType>&, Value>
{
   typedef void result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, RowMajor, InnerType, OuterType>& m,
                   size_type i, size_type j, Value const& x)
   {
      set_element(get_element(m.vec(),i), j, x);
   }
};

template <typename T, typename InnerType, typename OuterType, typename Value>
struct SetMatrixElement<SparseMatrix<T, ColMajor, InnerType, OuterType>&, Value>
{
   typedef void result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, ColMajor, InnerType, OuterType>& m,
                   size_type i, size_type j, Value const& x)
   {
      set_element(get_element(m.vec(),j), i, x);
   }
};

// get_element

template <typename SparseVec, typename Enable = void>
struct GetMatrixElement_Sparse {};

template <typename T, typename InnerType, typename OuterType>
struct GetMatrixElement_Sparse<
   SparseMatrix<T, RowMajor, InnerType, OuterType>
 , typename boost::enable_if<
      exists<
         typename GetVectorElement<
            typename reference_to_arg<
               typename GetVectorElement<OuterType>::result_type
            >::type
         >::result_type
      >
   >::type
>
{
   typedef typename GetVectorElement<OuterType>::result_type outer_result_type;
   typedef typename reference_to_arg<outer_result_type>::type outer_arg_type;
   typedef typename GetVectorElement<outer_arg_type>::result_type result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(get_element(m.vec(), i), j);
   }
};

template <typename T, typename InnerType, typename OuterType>
struct GetMatrixElement_Sparse<
   SparseMatrix<T, RowMajor, InnerType, OuterType>&
 , typename boost::enable_if<
      exists<
         typename GetVectorElement<
            typename GetVectorElement<typename boost::add_reference<OuterType>::type>::result_type
         >::result_type
      >
   >::type
>
{
   typedef typename GetVectorElement<typename boost::add_reference<OuterType>::type>::result_type
   outer_result_type;
   typedef typename GetVectorElement<outer_result_type>::result_type result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(get_element(m.vec(), i), j);
   }
};

template <typename T, typename InnerType, typename OuterType>
struct GetMatrixElement_Sparse<
   SparseMatrix<T, ColMajor, InnerType, OuterType>
 , typename boost::enable_if<
      exists<
         typename GetVectorElement<
            typename reference_to_arg<
               typename GetVectorElement<OuterType>::result_type
            >::type
         >::result_type
      >
   >::type
>
{
   typedef typename GetVectorElement<OuterType>::result_type outer_result_type;
   typedef typename reference_to_arg<outer_result_type>::type outer_arg_type;
   typedef typename GetVectorElement<outer_arg_type>::result_type result_type;
   typedef SparseMatrix<T, ColMajor, InnerType, OuterType> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(get_element(m.vec(), j), i);
   }
};

template <typename T, typename InnerType, typename OuterType>
struct GetMatrixElement_Sparse<
   SparseMatrix<T, ColMajor, InnerType, OuterType>&
 , typename boost::enable_if<
      exists<
         typename GetVectorElement<
            typename GetVectorElement<typename boost::add_reference<OuterType>::type>::result_type
         >::result_type
      >
   >::type
>
{
   typedef typename GetVectorElement<typename boost::add_reference<OuterType>::type>::result_type
   outer_result_type;
   typedef typename GetVectorElement<outer_result_type>::result_type result_type;
   typedef SparseMatrix<T, ColMajor, InnerType, OuterType>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(get_element(m.vec(), j), i);
   }
};

template <typename T, typename InnerType, typename OuterType, typename Orient>
struct GetMatrixElement<SparseMatrix<T, Orient, InnerType, OuterType> >
   : GetMatrixElement_Sparse<SparseMatrix<T, Orient, InnerType, OuterType> > {};

template <typename T, typename InnerType, typename OuterType, typename Orient>
struct GetMatrixElement<SparseMatrix<T, Orient, InnerType, OuterType>&>
   : GetMatrixElement_Sparse<SparseMatrix<T, Orient, InnerType, OuterType>&> {};

// add_element

template <typename T, typename InnerType, typename OuterType, typename Value>
struct AddMatrixElement<SparseMatrix<T, RowMajor, InnerType, OuterType>&, Value>
{
   typedef void result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, RowMajor, InnerType, OuterType>& m,
                   size_type i, size_type j, Value const& x) const
   {
      add_element(get_element(m.vec(),i), j, x);
   }
};

template <typename T, typename InnerType, typename OuterType, typename Value>
struct AddMatrixElement<SparseMatrix<T, ColMajor, InnerType, OuterType>&, Value>
{
   typedef void result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, ColMajor, InnerType, OuterType>& m,
                   size_type i, size_type j, Value const& x) const
   {
      add_element(get_element(m.vec(),j), i, x);
   }
};

template <typename T, typename InnerType, typename OuterType, typename Value, typename Float>
void add_element_cull(SparseMatrix<T, RowMajor, InnerType, OuterType>& m,
                      size_type i, size_type j, Value const& x,
                      Float const& Tol)
{
   add_element_cull(get_element(m.vec(),i), j, x, Tol);
}

// subtract_element

template <typename T, typename InnerType, typename OuterType, typename Value>
struct SubtractMatrixElement<SparseMatrix<T, RowMajor, InnerType, OuterType>&, Value>
{
   typedef void result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, ColMajor, InnerType, OuterType>& m,
                   size_type i, size_type j, Value const& x) const
   {
      subtract_element(get_element(m.vec(),i), j, x);
   }
};

template <typename T, typename InnerType, typename OuterType, typename Value>
struct SubtractMatrixElement<SparseMatrix<T, ColMajor, InnerType, OuterType>&, Value>
{
   typedef void result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, ColMajor, InnerType, OuterType>& m,
                   size_type i, size_type j, Value const& x) const
   {
      subtract_element(get_element(m.vec(),j), i, x);
   }
};

template <typename T, typename InnerType, typename OuterType, typename Value, typename Float>
void subtract_element_cull(SparseMatrix<T, RowMajor, InnerType, OuterType>& m,
                           size_type i, size_type j, Value const& x,
                           Float const& Tol)
{
   subtract_element_cull(get_element(m.vec(),j), i, x, Tol);
}

// zero_element

template <typename T, typename InnerType, typename OuterType>
struct ZeroMatrixElement<SparseMatrix<T, RowMajor, InnerType, OuterType>&>
{
   typedef void result_type;
   typedef SparseMatrix<T, RowMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, RowMajor, InnerType, OuterType>& m,
                   size_type i, size_type j) const
   {
      zero_element(get_element(m.vec(),i), j);
   }
};

template <typename T, typename InnerType, typename OuterType>
struct ZeroMatrixElement<SparseMatrix<T, ColMajor, InnerType, OuterType>&>
{
   typedef void result_type;
   typedef SparseMatrix<T, ColMajor, InnerType, OuterType>& first_argument_type;
   void operator()(SparseMatrix<T, ColMajor, InnerType, OuterType>& m,
                   size_type i, size_type j) const
   {
      zero_element(get_element(m.vec(),j), i);
   }
};

template <typename T, typename U, typename Orient, typename Nested>
struct MatrixDirectProduct<SparseMatrix<T, Orient>, SparseMatrix<U, Orient>, Nested>
{
   typedef typename make_value<typename Nested::result_type>::type result_value_type;
   typedef SparseMatrix<result_value_type, Orient> result_type;
   typedef SparseMatrix<T, Orient> const& first_argument_type;
   typedef SparseMatrix<U, Orient> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y, Nested f) const
   {
      size_type xs1 = size1(x);
      size_type xs2 = size2(x);
      size_type ys1 = size1(y);
      size_type ys2 = size2(y);
      typedef typename const_iterator<SparseMatrix<T, Orient> >::type x_iterator;
      typedef typename const_iterator<SparseMatrix<U, Orient> >::type y_iterator;
      typedef typename const_inner_iterator<SparseMatrix<T, Orient> >::type x_inner_iterator;
      typedef typename const_inner_iterator<SparseMatrix<U, Orient> >::type y_inner_iterator;

      result_type Result(xs1*ys1, xs2*ys2);

      for (x_iterator I = iterate(x); I; ++I)
      {
         for (x_inner_iterator J = iterate(I); J; ++J)
         {

            for (x_iterator K = iterate(y); K; ++K)
            {
               for (x_inner_iterator L = iterate(K); L; ++L)
               {

                  set_element(Result,
                              J.index1() * ys1 + L.index1(),
                              J.index2() * ys2 + L.index2(),
                              f(*J, *L));
               }
            }
         }
      }
      return Result;
   }
};

} // namespace LinearAlgebra

#endif
