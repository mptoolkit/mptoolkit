/* -*- C++ -*- $Id$

   Created 2005-03-03 Ian McCulloch

  TODO: only the RowMajor * ColMajor version is optimized.
  A RowMajor * RowMajor would be useful.

*/

#if !defined(COEFFICIENT_MULTIPLY_H_DJSHC84YT789YUFP89RHJT89PJHFP8943)
#define COEFFICIENT_MULTIPLY_H_DJSHC84YT789YUFP89RHJT89PJHFP8943

#include "coefficient_parallel_prod.h"
#include "matrix.h"
#include "sparsematrix.h"

namespace LinearAlgebra
{

template <typename CF>
struct BinderIJ
{
   typedef typename CF::result_type result_type;
   typedef size_type argument_type;

   BinderIJ(CF const& cf, size_type i, size_type j) : cf_(cf), i_(i), j_(j) {}

   result_type operator()(size_type k) const { return cf_(i_,k,j_); }

   CF const& cf_;
   size_type i_, j_;
};

template <typename CF>
inline
BinderIJ<CF> bindij(CF const& cf, size_type i, size_type j)
{
   return BinderIJ<CF>(cf, i, j);
}

//
// assign_coefficient_product2
//

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSInterface = typename interface<LHS>::type,
	  typename M1Interface = typename interface<M1>::type,
	  typename M2Interface = typename interface<M2>::type>
struct AssignCoefficient_Product2 {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested>
inline
void assign_coefficient_product2(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
{
   AssignCoefficient_Product2<LHS, M1, M2, CF, Nested>::apply(lhs, m1, m2, cf, f);
}

template <typename LHS, typename M1, typename M2,
          typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      MATRIX_EXPRESSION(LHSv, LHSi),
		      MATRIX_EXPRESSION(M1v, M1i),
		      MATRIX_EXPRESSION(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      assign_coefficient_product2(lhs, eval_expression(m1), eval_expression(m2), cf, f);
   }
};

// assign dense

template <typename LHS, typename M1, typename M2,
          typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, M1Orient, M1i),
		      DENSE_MATRIX(M2v, M2Orient, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_coefficient_product2(TempL, Temp1, Temp2, cf, f);
      assign(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2,
          typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      STRIDE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, M1Orient, M1i),
		      DENSE_MATRIX(M2v, M2Orient, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      assign_coefficient_product2(lhs, Temp1, Temp2, cf, f);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to AssignCoefficient_Product2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient, 
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename CF,
          typename Nested>
struct AssignCoefficient_Product2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2,
          typename CF,
          typename Nested>
struct AssignCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                         CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    set_element(lhs, I.index(), J.index(), 
                        coefficient_parallel_prod(*I, *J, bindij(cf, I.index(), J.index()), f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct AssignCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                         CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      assign_coefficient_product2(lhs, swap_sort_order(m1), m2, cf, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct AssignCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                         CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      assign_coefficient_product2(lhs, m1, swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct AssignCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                         CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested f)
   {
      assign_coefficient_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      STRIDE_MATRIX(LHSv, LHSOrient, LHSi),
		      STRIDE_MATRIX(M1v, M1Orient, M1i),
		      STRIDE_MATRIX(M2v, M2Orient, M2i)>
: AssignCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, 
                                    M1, M1Orient, 
                                    M2, M2Orient,
                                    CF, Nested> {};

// assign sparse

template <typename T>
T const& sick_show(T const& x)
{
   TRACE(typeid(x).name());
   return x;
}

template <typename T>
T sick_show(value_with_zero<T> const& x)
{
   TRACE(typeid(x).name());
   return x;
}

inline
int sick_show(value_with_zero<Matrix<std::complex<double> > > const& x)
{
   TRACE(typeid(x).name());
   TRACE(x.get().size1())(x.get().size2());
   CHECK(x.get().size1() != 0);
   return 1;
}

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename CF, typename Nested>
struct AssignCoefficient_Product2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AssignCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      zero_all(lhs);
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element_check_if_zero(lhs, I.index(), J.index(), 
                   coefficient_parallel_prod(*I, *J, 
                                             bindij(cf, I.index(), J.index()), f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AssignCoefficient_Product2_Sparse<LHS, M1, ColMajor, M2, ColMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_coefficient_product2(lhs, Temp1, m2, cf, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AssignCoefficient_Product2_Sparse<LHS, M1, ColMajor, M2, RowMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_product2(lhs, Temp1, Temp2, cf, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AssignCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, RowMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_product2(lhs, m1, Temp2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      LOCAL_MATRIX(LHSv, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, M1Orient, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, M2Orient, M2i)>
: AssignCoefficient_Product2_Sparse<LHS, M1, M1Orient, M2, M2Orient, CF, Nested> {};

// assign mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, RowMajor, M1i),
		      DENSE_MATRIX(M2v, ColMajor, M2i)>
: AssignCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, ColMajor, M1i),
		      DENSE_MATRIX(M2v, ColMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_coefficient_product2(lhs, Temp1, m2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, ColMajor, M1i),
		      DENSE_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      assign_coefficient_product2(lhs, Temp1, swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, RowMajor, M1i),
		      DENSE_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      assign_coefficient_product2(lhs, m1, swap_sort_order(m2), cf, f);
   }
};


template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, RowMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, ColMajor, M2i)>
: AssignCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, ColMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, ColMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      assign_coefficient_product2(lhs, swap_sort_order(m1), m2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, ColMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_product2(lhs, swap_sort_order(m1), Temp2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, RowMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_coefficient_product2(lhs, m1, Temp2, cf, f);
   }
};

// assign mixed sparse/diagonal

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
				  SPARSE_MATRIX(LHSv, LHSi),
				  SPARSE_MATRIX(M1v, M1i),
				  DIAGONAL_MATRIX(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      zero_all(lhs);
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_inner_iterator<M1>::type J = iterate(I);
	 while (J)
	 {
	    add_element(lhs, J.index1(), J.index2(), 
			cf(J.index1(), J.index2(), J.index2()) * f(*J, m2.diagonal()[J.index2()]));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AssignCoefficient_Product2<LHS, M1, M2, CF, Nested,
				  SPARSE_MATRIX(LHSv, LHSi),
				  DIAGONAL_MATRIX(M1v, M1i),
				  SPARSE_MATRIX(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      zero_all(lhs);
      typename const_iterator<M2>::type I = iterate(m2);
      while (I)
      {
	 typename const_inner_iterator<M2>::type J = iterate(I);
	 while (J)
	 {
	    add_element(lhs, J.index1(), J.index2(), 
			cf(J.index1(), J.index1(), J.index2()) * f(m1.diagonal()[J.index1()], *J));
	    ++J;
	 }
	 ++I;
      }
   }
};

//
// add_coefficient_product2
//

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSInterface = typename interface<LHS>::type,
	  typename M1Interface = typename interface<M1>::type,
	  typename M2Interface = typename interface<M2>::type>
struct AddCoefficient_Product2 {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested>
inline
void add_coefficient_product2(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
{
   AddCoefficient_Product2<LHS, M1, M2, CF, Nested>::apply(lhs, m1, m2, cf, f);
}

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      MATRIX_EXPRESSION(LHSv, LHSi),
		      MATRIX_EXPRESSION(M1v, M1i),
		      MATRIX_EXPRESSION(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      add_coefficient_product2(lhs, eval_expression(m1), eval_expression(m2), cf, f);
   }
};

// add dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, M1Orient, M1i),
		      DENSE_MATRIX(M2v, M2Orient, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_coefficient_product2(TempL, Temp1, Temp2, cf, f);
      add(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      STRIDE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, M1Orient, M1i),
		      DENSE_MATRIX(M2v, M2Orient, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      add_coefficient_product2(lhs, Temp1, Temp2, cf, f);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to AddCoefficient_Product2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient, 
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename CF, typename Nested>
struct AddCoefficient_Product2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2,
          typename CF, typename Nested>
struct AddCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element(lhs, I.index(), J.index(), 
                        coefficient_parallel_prod(*I, *J, bindij(cf, I.index(), J.index()), f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct AddCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      add_coefficient_product2(lhs, swap_sort_order(m1), m2, cf, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct AddCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      add_coefficient_product2(lhs, m1, swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct AddCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      add_coefficient_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      STRIDE_MATRIX(LHSv, LHSOrient, LHSi),
		      STRIDE_MATRIX(M1v, M1Orient, M1i),
		      STRIDE_MATRIX(M2v, M2Orient, M2i)>
: AddCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, M1Orient, M2, M2Orient,
                                             CF, Nested> {};


// add sparse

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient, typename CF, typename Nested>
struct AddCoefficient_Product2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AddCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element_check_if_zero(lhs, I.index(), J.index(), 
                 coefficient_parallel_prod(*I, *J, bindij(cf, I.index(), J.index()), f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AddCoefficient_Product2_Sparse<LHS, M1, ColMajor, M2, ColMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_coefficient_product2(lhs, Temp1, m2, cf, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AddCoefficient_Product2_Sparse<LHS, M1, ColMajor, M2, RowMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_product2(lhs, Temp1, Temp2, cf, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct AddCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, RowMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_product2(lhs, m1, Temp2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      LOCAL_MATRIX(LHSv, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, M1Orient, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, M2Orient, M2i)>
: AddCoefficient_Product2_Sparse<LHS, M1, M1Orient, M2, M2Orient, CF, Nested> {};

// add mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, RowMajor, M1i),
		      DENSE_MATRIX(M2v, ColMajor, M2i)>
: AddCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, ColMajor, M1i),
		      DENSE_MATRIX(M2v, ColMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_coefficient_product2(lhs, Temp1, m2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, ColMajor, M1i),
		      DENSE_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_coefficient_product2(lhs, Temp1, swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, RowMajor, M1i),
		      DENSE_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      add_coefficient_product2(lhs, m1, swap_sort_order(m2), cf, f);
   }
};


template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, RowMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, ColMajor, M2i)>
: AddCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, ColMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, ColMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      add_coefficient_product2(lhs, swap_sort_order(m1), m2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, ColMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_product2(lhs, swap_sort_order(m1), Temp2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, RowMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_coefficient_product2(lhs, m1, Temp2, cf, f);
   }
};

// add mixed sparse/diagonal

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
			       SPARSE_MATRIX(LHSv, LHSi),
			       SPARSE_MATRIX(M1v, M1i),
			       DIAGONAL_MATRIX(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_inner_iterator<M1>::type J = iterate(I);
	 while (J)
	 {
	    add_element(lhs, J.index1(), J.index2(), 
			cf(J.index1(), J.index2(), J.index2()) * f(*J, m2.diagonal()[J.index2()]));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct AddCoefficient_Product2<LHS, M1, M2, CF, Nested,
			       SPARSE_MATRIX(LHSv, LHSi),
			       DIAGONAL_MATRIX(M1v, M1i),
			       SPARSE_MATRIX(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M2>::type I = iterate(m2);
      while (I)
      {
	 typename const_inner_iterator<M2>::type J = iterate(I);
	 while (J)
	 {
	    add_element(lhs, J.index1(), J.index2(), 
			cf(J.index1(), J.index1(), J.index2()) * f(m1.diagonal()[J.index1()], *J));
	    ++J;
	 }
	 ++I;
      }
   }
};

//
// subtract_coefficient_product2
//

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSInterface = typename interface<LHS>::type,
	  typename M1Interface = typename interface<M1>::type,
	  typename M2Interface = typename interface<M2>::type>
struct SubtractCoefficient_Product2 {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested>
inline
void subtract_coefficient_product2(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
{
   SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested>::apply(lhs, m1, m2, cf, f);
}

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      MATRIX_EXPRESSION(LHSv, LHSi),
		      MATRIX_EXPRESSION(M1v, M1i),
		      MATRIX_EXPRESSION(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      subtract_coefficient_product2(lhs, eval_expression(m1), eval_expression(m2), cf, f);
   }
};

// subtract dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, M1Orient, M1i),
		      DENSE_MATRIX(M2v, M2Orient, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      Matrix<LHSv, RowMajor> TempL(size1(m1), size2(m2));
      assign_coefficient_product2(TempL, Temp1, Temp2, cf, f);
      subtract(lhs, TempL);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      STRIDE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, M1Orient, M1i),
		      DENSE_MATRIX(M2v, M2Orient, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      Matrix<M1v, RowMajor> Temp1(m1);
      Matrix<M2v, ColMajor> Temp2(m2);
      subtract_coefficient_product2(lhs, Temp1, Temp2, cf, f);
   }
};

// a partial specialization for LHS = RowMajor * ColMajor
// would be rude here - it would force specializations on the
// LHS or the value_type to also follow suit.
// Instead, we forward to SubtractCoefficient_Product2_StrideStrideStride
// and specialize further there.

template <typename LHS, typename LHSOrient, 
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename CF, typename Nested>
struct SubtractCoefficient_Product2_StrideStrideStride;

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2,
          typename CF, typename Nested>
struct SubtractCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, ColMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    subtract_element(lhs, I.index(), J.index(), 
                        coefficient_parallel_prod(*I, *J, bindij(cf, I.index(), J.index()), f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct SubtractCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, ColMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      subtract_coefficient_product2(lhs, swap_sort_order(m1), m2, cf, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct SubtractCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, RowMajor, M2, RowMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      subtract_coefficient_product2(lhs, m1, swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename LHSOrient,
	  typename M1, typename M2, typename CF, typename Nested>
struct SubtractCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, ColMajor, M2, RowMajor,
                                      CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      subtract_coefficient_product2(lhs, swap_sort_order(m1), swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      STRIDE_MATRIX(LHSv, LHSOrient, LHSi),
		      STRIDE_MATRIX(M1v, M1Orient, M1i),
		      STRIDE_MATRIX(M2v, M2Orient, M2i)>
: SubtractCoefficient_Product2_StrideStrideStride<LHS, LHSOrient, M1, M1Orient, M2, M2Orient,
                                             CF, Nested> {};

// subtract sparse

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient, typename CF, typename Nested>
struct SubtractCoefficient_Product2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct SubtractCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    subtract_element_check_if_zero(lhs, I.index(), J.index(), 
                 coefficient_parallel_prod(*I, *J, bindij(cf, I.index(), J.index()), f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct SubtractCoefficient_Product2_Sparse<LHS, M1, ColMajor, M2, ColMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_coefficient_product2(lhs, Temp1, m2, cf, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct SubtractCoefficient_Product2_Sparse<LHS, M1, ColMajor, M2, RowMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_product2(lhs, Temp1, Temp2, cf, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename CF, typename Nested>
struct SubtractCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, RowMajor, CF, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_product2(lhs, m1, Temp2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      LOCAL_MATRIX(LHSv, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, M1Orient, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, M2Orient, M2i)>
: SubtractCoefficient_Product2_Sparse<LHS, M1, M1Orient, M2, M2Orient, CF, Nested> {};


// subtract mixed sparse/dense

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, RowMajor, M1i),
		      DENSE_MATRIX(M2v, ColMajor, M2i)>
: SubtractCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, ColMajor, M1i),
		      DENSE_MATRIX(M2v, ColMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_coefficient_product2(lhs, Temp1, m2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, ColMajor, M1i),
		      DENSE_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      subtract_coefficient_product2(lhs, Temp1, swap_sort_order(m2), cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, RowMajor, M1i),
		      DENSE_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      subtract_coefficient_product2(lhs, m1, swap_sort_order(m2), cf, f);
   }
};


template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, RowMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, ColMajor, M2i)>
: SubtractCoefficient_Product2_Sparse<LHS, M1, RowMajor, M2, ColMajor, CF, Nested> {};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, ColMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, ColMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      subtract_coefficient_product2(lhs, swap_sort_order(m1), m2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, ColMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_product2(lhs, swap_sort_order(m1), Temp2, cf, f);
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSOrient, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
		      DENSE_MATRIX(LHSv, LHSOrient, LHSi),
		      DENSE_MATRIX(M1v, RowMajor, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, RowMajor, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      subtract_coefficient_product2(lhs, m1, Temp2, cf, f);
   }
};

// subtract mixed sparse/diagonal

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
			       SPARSE_MATRIX(LHSv, LHSi),
			       SPARSE_MATRIX(M1v, M1i),
			       DIAGONAL_MATRIX(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_inner_iterator<M1>::type J = iterate(I);
	 while (J)
	 {
	    subtract_element(lhs, J.index1(), J.index2(), 
			     cf(J.index1(), J.index2(), J.index2()) * f(*J, m2.diagonal()[J.index2()]));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS, typename M1, typename M2, typename CF, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1i,
	  typename M2v, typename M2i>
struct SubtractCoefficient_Product2<LHS, M1, M2, CF, Nested,
			       SPARSE_MATRIX(LHSv, LHSi),
			       DIAGONAL_MATRIX(M1v, M1i),
			       SPARSE_MATRIX(M2v, M2i)>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, CF const& cf, Nested const& f)
   {
      typename const_iterator<M2>::type I = iterate(m2);
      while (I)
      {
	 typename const_inner_iterator<M2>::type J = iterate(I);
	 while (J)
	 {
	    subtract_element(lhs, J.index1(), J.index2(), 
			     cf(J.index1(), J.index1(), J.index2()) * f(m1.diagonal()[J.index1()], *J));
	    ++J;
	 }
	 ++I;
      }
   }
};

} // namespace LinearAlgebra

#endif
