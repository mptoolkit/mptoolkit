/* -*- C++ -*- $Id$

  sparse matrix*matrix operations

  Created 2005-02-23 Ian McCulloch
*/

#if !defined(SPARSEMATRIXPRODUCT_H_DJSHC84YT789YUFP89RHJT89PJHFP8943)
#define SPARSEMATRIXPRODUCT_H_DJSHC84YT789YUFP89RHJT89PJHFP8943

#include "sparsematrix.h"
#include "vector.h"
#include "mapvector.h"

namespace LinearAlgebra
{

//
// assign_product2
//

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient,
          typename Nested>
struct AssignProduct2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      zero_all(lhs);
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element_check_if_zero(lhs, I.index(), J.index(), parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, ColMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp(m1);
      assign_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, ColMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_product2(lhs, Temp1, Temp2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AssignProduct2_Sparse<LHS, M1, RowMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      assign_product2(lhs, m1, Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Nested,
		      LOCAL_MATRIX(LHSv, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, M1Orient, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, M2Orient, M2i)>
: AssignProduct2_Sparse<LHS, M1, M1Orient, M2, M2Orient, Nested> {};

//
// add_product2
//

template <typename LHS,
	  typename M1, typename M1Orient,
	  typename M2, typename M2Orient, typename Nested>
struct AddProduct2_Sparse;

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, RowMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      typename const_iterator<M1>::type I = iterate(m1);
      while (I)
      {
	 typename const_iterator<M2>::type J = iterate(m2);
	 while (J)
	 {
	    add_element_check_if_zero(lhs, I.index(), J.index(), parallel_prod(*I, *J, f));
	    ++J;
	 }
	 ++I;
      }
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, ColMajor, M2, ColMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      add_product2(lhs, Temp1, m2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, ColMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M1>::value_type, RowMajor> Temp1(m1);
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_product2(lhs, Temp1, Temp2, f);
   }
};

template <typename LHS,
	  typename M1, typename M2, typename Nested>
struct AddProduct2_Sparse<LHS, M1, RowMajor, M2, RowMajor, Nested>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, Nested const& f)
   {
      SparseMatrix<typename interface<M2>::value_type, ColMajor> Temp2(m2);
      add_product2(lhs, m1, Temp2, f);
   }
};

template <typename LHS, typename M1, typename M2, typename Nested,
	  typename LHSv, typename LHSi,
	  typename M1v, typename M1Orient, typename M1i,
	  typename M2v, typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Nested,
		      LOCAL_MATRIX(LHSv, LHSi),
		      COMPRESSED_OUTER_MATRIX(M1v, M1Orient, M1i),
		      COMPRESSED_OUTER_MATRIX(M2v, M2Orient, M2i)>
: AddProduct2_Sparse<LHS, M1, M1Orient, M2, M2Orient, Nested> {};

} // namespace LinearAlgebra

#endif
