/* -*- C++ -*- $Id$

  compressedmatrix.h

  sparsematrix.h makes most of this header redundant.
  The only part of this header that is still in use is the CompressedMatrixBase
  helper class.  TODO: something like the CompressedMatrix interface is needed to
  construct a 'matrix view' of a vector of vectors.

*/

#if !defined(COMPRESSEDMATRIX_H_A636DJBNDJ48976BNM6NM64RJYG86V5)
#define COMPRESSEDMATRIX_H_A636DJBNDJ48976BNM6NM64RJYG86V5

#include "matrixinterface.h"
#include "matrixoperationsbase.h"
#include "mapvector.h"        // to construct the 'default' sparse type
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif
#include "matrixiterators.h"
#include "crtp_matrix.h"
#include <boost/mpl/assert.hpp>

namespace LinearAlgebra
{

template <typename VecOfVec, typename Derived>
class CompressedMatrixBase : public MatrixBase<Derived>
{
   public:
      BOOST_MPL_ASSERT((is_vector<VecOfVec>));

      typedef VecOfVec data_type;
      typedef typename interface<VecOfVec>::type outer_interface;
      typedef typename interface<VecOfVec>::value_type outer_value_type;

      BOOST_MPL_ASSERT((is_vector<outer_value_type>));

      typedef typename interface<outer_value_type>::type inner_interface;
      typedef typename interface<outer_value_type>::value_type value_type;

      typedef typename GetVectorElement<data_type>::result_type  const_outer_reference;

   //      typedef typename GetElement<typename 
   //         basic_type<outer_reference>::type&>::result_type reference;

      typedef typename GetVectorElement<typename 
         basic_type<const_outer_reference>::type>::result_type const_reference;

};

template <typename VecOfVec, typename Orientation = RowMajor>
class CompressedMatrix;

template <typename VecOfVec>
class CompressedMatrix<VecOfVec, RowMajor> 
   : public CompressedMatrixBase<VecOfVec, CompressedMatrix<VecOfVec, RowMajor> >
{
   public:
      typedef CompressedMatrixBase<VecOfVec, CompressedMatrix<VecOfVec, RowMajor> > base_type;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::outer_value_type outer_value_type;

      typedef typename make_reference<VecOfVec>::type       vec_reference;
      typedef typename make_const_reference<VecOfVec>::type const_vec_reference;

      CompressedMatrix(size_t Rows, size_t Cols)
	 : InnerSize_(Cols), Data_(Rows, outer_value_type(Cols)) {}

      CompressedMatrix(CompressedMatrix const& x) 
	 : InnerSize_(x.inner_size()), Data_(x.Data_) {}

      template <typename U>
      CompressedMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type* dummy = 0);

      template <typename Arg1, typename Arg2>
      typename MatrixBracket<CompressedMatrix, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y) const
      {
	 return MatrixBracket<CompressedMatrix, Arg1, Arg2>()(*this, x,y);
      }

      template <typename Arg1, typename Arg2>
      typename MatrixBracket<CompressedMatrix&, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y)
      {
	 return MatrixBracket<CompressedMatrix&, Arg1, Arg2>()(*this, x,y);
      }

      vec_reference vec() { return Data_; }
      const_vec_reference vec() const { return Data_; }

      size_type inner_size() const { return InnerSize_; }

   private:
      size_type InnerSize_;
      VecOfVec Data_;
};

template <typename VecOfVec>
class CompressedMatrix<VecOfVec, ColMajor> 
   : public CompressedMatrixBase<VecOfVec, CompressedMatrix<VecOfVec, ColMajor> >
{
   public:
      typedef CompressedMatrixBase<VecOfVec, CompressedMatrix<VecOfVec, ColMajor> > base_type;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::outer_value_type outer_value_type;

      typedef typename make_reference<VecOfVec>::type       vec_reference;
      typedef typename make_const_reference<VecOfVec>::type const_vec_reference;

      CompressedMatrix(size_t Rows, size_t Cols)
	 : InnerSize_(Rows), Data_(Cols, outer_value_type(Rows)) {}

      CompressedMatrix(CompressedMatrix const& x) 
	 : InnerSize_(x.inner_size()), Data_(x.Data_) {}

      template <typename U>
      CompressedMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type* dummy = 0);

      template <typename Arg1, typename Arg2>
      typename MatrixBracket<CompressedMatrix, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y) const
      {
	 return MatrixBracket<CompressedMatrix, Arg1, Arg2>()(*this, x,y);
      }

      template <typename Arg1, typename Arg2>
      typename MatrixBracket<CompressedMatrix&, Arg1, Arg2>::result_type
      operator()(Arg1 const& x, Arg2 const& y)
      {
	 return MatrixBracket<CompressedMatrix&, Arg1, Arg2>()(*this, x,y);
      }

      vec_reference vec() { return Data_; }
      const_vec_reference vec() const { return Data_; }

      size_type inner_size() const { return InnerSize_; }

   private:
      size_type InnerSize_;
      VecOfVec Data_;
};

template <typename VecOfVec>
template <typename U>
CompressedMatrix<VecOfVec, RowMajor>::
CompressedMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type*)
   : InnerSize_(size2(x)), Data_(size1(x), outer_value_type(InnerSize_))
{
   assign(*this, x);
}

template <typename VecOfVec>
template <typename U>
CompressedMatrix<VecOfVec, ColMajor>::
CompressedMatrix(U const& x, typename boost::enable_if<is_matrix<U> >::type*)
   : InnerSize_(size1(x)), Data_(size2(x), outer_value_type(InnerSize_))
{
   assign(*this, x);
}

// interface

template <typename VecOfVec, typename Orientation>
struct interface<CompressedMatrix<VecOfVec, Orientation> >
{
   typedef typename CompressedMatrix<VecOfVec, Orientation>::value_type value_type;
   typedef typename CompressedMatrix<VecOfVec, Orientation>::outer_interface Outer_;

   typedef COMPRESSED_OUTER_MATRIX_V(value_type, Orientation, Outer_, void) type;
};

// iterators

template <typename VecOfVec, typename Orientation>
struct Iterate<CompressedMatrix<VecOfVec, Orientation>&>
{
   typedef typename Iterate<VecOfVec&>::result_type iter;
   typedef MatrixOuterIterator<iter, Orientation> result_type;
   typedef CompressedMatrix<VecOfVec, Orientation>& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.vec()));
   }
};

template <typename VecOfVec, typename Orientation>
struct Iterate<CompressedMatrix<VecOfVec, Orientation> >
{
   typedef typename Iterate<VecOfVec>::result_type iter;
   typedef MatrixOuterIterator<iter, Orientation> result_type;
   typedef CompressedMatrix<VecOfVec, Orientation> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.vec()));
   }
};

// size

template <typename VecOfVec>
struct Size1<CompressedMatrix<VecOfVec, RowMajor> >
{
   typedef size_type result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor> const& argument_type;
   result_type operator()(argument_type x) const { return size(x.vec()); }
};

template <typename VecOfVec>
struct Size2<CompressedMatrix<VecOfVec, RowMajor> >
{
   typedef size_type result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor> const& argument_type;
   result_type operator()(argument_type x) const 
   { 
      return x.inner_size();
   }
};

template <typename VecOfVec>
struct Size2<CompressedMatrix<VecOfVec, ColMajor> >
{
   typedef size_type result_type;
   typedef CompressedMatrix<VecOfVec, ColMajor> const& argument_type;
   result_type operator()(argument_type x) const { return size(x.vec()); }
};

template <typename VecOfVec>
struct Size1<CompressedMatrix<VecOfVec, ColMajor> >
{
   typedef size_type result_type;
   typedef CompressedMatrix<VecOfVec, ColMajor> const& argument_type;
   result_type operator()(argument_type x) const 
   { 
      return x.inner_size();
   }
};

// vector_view

template <typename VecOfVec, typename Orient>
struct VectorView<CompressedMatrix<VecOfVec, Orient> >
{
   typedef typename make_const_reference<VecOfVec>::type result_type;
   typedef CompressedMatrix<VecOfVec, Orient> const& argument_type;
   result_type operator()(argument_type x) const { return x.vec(); }
};

template <typename VecOfVec, typename Orient>
struct VectorView<CompressedMatrix<VecOfVec, Orient>&>
{
   typedef typename make_reference<VecOfVec>::type result_type;
   typedef CompressedMatrix<VecOfVec, Orient>& argument_type;
   result_type operator()(argument_type x) const { return x.vec(); }
};

// get_element

template <typename CompressedVec, typename Enable = void>
struct GetMatrixElement_Compressed {};

template <typename VecOfVec>
struct GetMatrixElement_Compressed<
   CompressedMatrix<VecOfVec, RowMajor>
 , typename boost::enable_if<
      exists<
         typename GetVectorElement<
            typename reference_to_arg<
               typename GetVectorElement<VecOfVec>::result_type
            >::type
         >::result_type
      >
   >::type
>
{
   typedef typename GetVectorElement<VecOfVec>::result_type outer_result_type;
   typedef typename reference_to_arg<outer_result_type>::type outer_arg_type;
   typedef typename GetVectorElement<outer_arg_type>::result_type result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(get_element(m.vec(), i), j);
   }
};

template <typename VecOfVec>
struct GetMatrixElement_Compressed<
   CompressedMatrix<VecOfVec, RowMajor>&
 , typename boost::enable_if<
      exists<
         typename GetVectorElement<
            typename GetVectorElement<typename boost::add_reference<VecOfVec>::type>::result_type
         >::result_type
      >
   >::type
>
{
   typedef typename GetVectorElement<typename boost::add_reference<VecOfVec>::type>::result_type 
   outer_result_type;
   typedef typename GetVectorElement<outer_result_type>::result_type result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;

   result_type operator()(first_argument_type m, size_type i, size_type j) const
   {
      return get_element(get_element(m.vec(), i), j);
   }
};

template <typename VecOfVec, typename Orient>
struct GetMatrixElement<CompressedMatrix<VecOfVec, Orient> > 
   : GetMatrixElement_Compressed<CompressedMatrix<VecOfVec, Orient> > {};

template <typename VecOfVec, typename Orient>
struct GetMatrixElement<CompressedMatrix<VecOfVec, Orient>&> 
   : GetMatrixElement_Compressed<CompressedMatrix<VecOfVec, Orient>&> {};

// add_element

template <typename VecOfVec, typename Value>
struct AddMatrixElement<CompressedMatrix<VecOfVec, RowMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, RowMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      add_element(get_element(m.vec(),i), j, x);
   }
};

template <typename VecOfVec, typename Value>
struct AddMatrixElement<CompressedMatrix<VecOfVec, ColMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, ColMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      add_element(get_element(m.vec(),j), i, x);
   }
};

// subtract_element

template <typename VecOfVec, typename Value>
struct SubtractMatrixElement<CompressedMatrix<VecOfVec, RowMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, ColMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      subtract_element(get_element(m.vec(),i), j, x);
   }
};

template <typename VecOfVec, typename Value>
struct SubtractMatrixElement<CompressedMatrix<VecOfVec, ColMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, ColMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      subtract_element(get_element(m.vec(),j), i, x);
   }
};

// zero_element

template <typename VecOfVec, typename Value>
struct ZeroMatrixElement<CompressedMatrix<VecOfVec, RowMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, ColMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      zero_element(get_element(m.vec(),i), j, x);
   }
};

template <typename VecOfVec, typename Value>
struct ZeroMatrixElement<CompressedMatrix<VecOfVec, ColMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, ColMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      zero_element(get_element(m.vec(),j), i, x);
   }
};

// set_element

template <typename VecOfVec, typename Value>
struct SetMatrixElement<CompressedMatrix<VecOfVec, RowMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, ColMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      set_element(get_element(m.vec(),i), j, x);
   }
};

template <typename VecOfVec, typename Value>
struct SetMatrixElement<CompressedMatrix<VecOfVec, ColMajor>, Value> 
{
   typedef void result_type;
   typedef CompressedMatrix<VecOfVec, RowMajor>& first_argument_type;
   void operator()(CompressedMatrix<VecOfVec, ColMajor>& m, 
		   size_type i, size_type j, Value const& x)
   {
      set_element(get_element(m.vec(),j), i, x);
   }
};

// zero_all

template <typename VecOfVec, typename Orient>
struct ZeroAll<CompressedMatrix<VecOfVec, Orient> >
{
   static void apply( CompressedMatrix<VecOfVec, Orient>& m)
   {
      return zero_all(m.vec());
   }
};

} // namespace LinearAlgebra

#endif
