/* -*- C++ -*- $Id$
   
  matrixmatrixmultiplication.h 

  Created 2005-02-21 Ian McCulloch
*/

#if !defined(MATRIXMATRIXMULTIPLICATION_H_HJUIGH579879RHEOP)
#define MATRIXMATRIXMULTIPLICATION_H_HJUIGH579879RHEOP

#include "matrixproductoperations.h"
#include <boost/utility/result_of.hpp>

namespace LinearAlgebra
{

//
// MatrixProductProxy
//

template <typename T1, typename T2, typename NestedMult>
class MatrixProductProxy;

template <typename T1, typename T2, typename NestedMult>
struct abstract_interface<MatrixProductProxy<T1, T2, NestedMult> >
{
   typedef typename matrix_abstract_or<abstract_interface<T1>,
                                       abstract_interface<T2> >::type type;
};

template <typename T1, typename T2, typename NestedMult>
class MatrixProductProxy
{
   public:
      BOOST_MPL_ASSERT((is_const_proxy_reference<T1>));
      BOOST_MPL_ASSERT((is_const_proxy_reference<T2>));

      // mark this type as a const proxy reference
      typedef boost::mpl::true_ const_proxy;

      typedef T1 reference1;
      typedef T2 reference2;

      typedef typename basic_type<T1>::type matrix1_type;
      typedef typename basic_type<T2>::type matrix2_type;

      typedef NestedMult functor_type;

      size_type size1() const { return Size1<matrix1_type>()(x_); }
      size_type size2() const { return Size2<matrix2_type>()(y_); }

      reference1 matrix1() const { return x_; }
      reference2 matrix2() const { return y_; }

      MatrixProductProxy(reference1 x, reference2 y, 
                         functor_type f = functor_type()) : x_(x), y_(y), f_(f)
      {
	 using LinearAlgebra::size1; using LinearAlgebra::size2;
	 DEBUG_PRECONDITION_EQUAL(size2(x), size1(y));
      }

      functor_type functor() const { return f_; }

   private:
      MatrixProductProxy& operator=(MatrixProductProxy const&); // not implemented

      reference1 x_;
      reference2 y_;

      functor_type f_;
};

// interface

template <typename T1, typename T2, typename Nested>
struct interface<MatrixProductProxy<T1, T2, Nested> > 
{
   typedef typename interface<T1>::value_type value1_type;
   typedef typename interface<T2>::value_type value2_type;
   //   typedef typename make_value<typename Nested::result_type>::type value_type;
   typedef typename make_value<
      typename boost::result_of<
         Nested(value1_type, value2_type)
      >::type
   >::type value_type;

   typedef Concepts::MatrixExpression<value_type, void> type;
};

// Assign

template <typename LHS, typename T1, typename T2, typename F>
struct AssignExpression<LHS&, MatrixProductProxy<T1, T2, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixProductProxy<T1, T2, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      assign_product2(x, y.matrix1(), y.matrix2(), y.functor());
   }
};

// Add

template <typename LHS, typename T1, typename T2, typename F>
struct AddExpression<LHS&, MatrixProductProxy<T1, T2, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixProductProxy<T1, T2, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      add_product2(x, y.matrix1(), y.matrix2(), y.functor());
   }
};

// Subtract

template <typename LHS, typename T1, typename T2, typename F>
struct SubtractExpression<LHS&, MatrixProductProxy<T1, T2, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixProductProxy<T1, T2, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      subtract_product2(x, y.matrix1(), y.matrix2(), y.functor());
   }
};

// 
// MatrixMatrixMultiplication
//

template <typename M1, typename M2, 
          typename Nested,
          typename M1v, typename M1i, 
          typename M2v, typename M2i>
struct MatrixMatrixMultiplication<M1, M2, Nested,
				  Concepts::AnyMatrix<M1v, M1i>, 
				  Concepts::AnyMatrix<M2v, M2i>>
{
   typedef M1 const& first_argument_type;
   typedef M2 const& second_argument_type;
   typedef MatrixProductProxy<typename make_const_reference<M1>::type,
			      typename make_const_reference<M2>::type,
                              Nested> result_type;

   result_type operator()(M1 const& m1, M2 const& m2) const { return result_type(m1, m2); }

   result_type operator()(M1 const& m1, M2 const& m2, Nested const& n) const { return result_type(m1, m2, n); }
};

// transpose

template <typename M1, typename M2>
struct Transpose<MatrixProductProxy<M1, M2, 
                                    Multiplication<typename basic_type<M1>::type::value_type,
                                                   typename basic_type<M2>::type::value_type> > >
{
   typedef typename basic_type<M1>::type M1Type;
   typedef typename basic_type<M2>::type M2Type;
   typedef typename Transpose<M1Type>::result_type TransM1Type;
   typedef typename Transpose<M2Type>::result_type TransM2Type;
   typedef MatrixProductProxy<TransM2Type, TransM1Type, Multiplication<typename M2Type::value_type,
                                                                       typename M1Type::value_type> >
   result_type;

   typedef typename boost::mpl::print<M1>::type dummy;

   typedef MatrixProductProxy<M1, M2, Multiplication<typename basic_type<M1>::type::value_type,
                                                   typename basic_type<M2>::type::value_type> >
   argument_type;

   result_type operator()(argument_type x) const
   {
      return result_type(trans(x.matrix2()), trans(x.matrix1()));
   }
};

} // namespace LinearAlgebra

#endif
