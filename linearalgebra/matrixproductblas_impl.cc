/* -*- C++ -*- $Id$

  matrixproductblas_impl.cc

  Created 2005-02-21 Ian McCulloch
*/

#if !defined(MATRIXPRODUCTBLAS_IMPL_H_DSHJFW5Y789FY78943FOH)
#define MATRIXPRODUCTBLAS_IMPL_H_DSHJFW5Y789FY78943FOH

#include "common/blas3f.h"
#include "matrixproductoperations.h"

namespace LinearAlgebra
{

//
// real assign
//

template <typename LHS, typename M1, typename M2,
	  typename LHSOrient, typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::StrideMatrix<double, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<double, double> f);
};

template <typename LHS, typename M1, typename M2,
	  typename LHSOrient, typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
void
AssignProduct2<LHS, M1, M2, Multiplication<double, double>,
               Concepts::StrideMatrix<double, LHSOrient, LHSi>,
               Concepts::StrideMatrix<double, M1Orient, M1i>,
               Concepts::StrideMatrix<double, M2Orient, M2i>>
::apply(LHS& lhs, M1 const& m1, M2 const& m2,
        Multiplication<double, double> f)
{
   if (!is_blas_matrix(lhs))
   {
      Matrix<double, RowMajor> TempL(size1(m1), size2(m2));
      if (!is_blas_matrix(m1))
      {
         Matrix<double, RowMajor> Temp1(m1);
         if (is_blas_matrix(m2))
	 {
            assign_product2(TempL, Temp1, m2);
         }
         else
	 {
            Matrix<double, ColMajor> Temp2(m2);
            assign_product2(TempL, Temp1, Temp2);
         }
         assign(lhs, TempL);
      }
      return;
   }
   // else
   if (!is_blas_matrix(m1))
   {
      Matrix<double, RowMajor> Temp1(m1);
      if (is_blas_matrix(m2))
      {
         assign_product2(lhs, Temp1, m2);
      }
      else
      {
         Matrix<double, ColMajor> Temp2(m2);
         assign_product2(lhs, Temp1, Temp2);
      }
      return;
   }
   // else
   if (!is_blas_matrix(m2))
   {
      Matrix<double, ColMajor> Temp2(m2);
      assign_product2(lhs, m1, Temp2);
      return;
   }
   // else
   DEBUG_CHECK(is_blas_matrix(lhs) && is_blas_matrix(m1) && is_blas_matrix(m2));
   if (blas_trans_col(lhs) == 'T')
   {
      BLAS::dgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
                  1, data(m2), leading_dimension(m2), 
                  data(m1), leading_dimension(m1),
                  0, data(lhs), leading_dimension(lhs));
   }
   else
   {
      BLAS::dgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
                  1, data(m1), leading_dimension(m1), 
                  data(m2), leading_dimension(m2),
                  0, data(lhs), leading_dimension(lhs));
   }
}

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::ContiguousMatrix<double, RowMajor, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<double, double> f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      //      DEBUG_CHECK_EQUAL(stride2(lhs), 1);
      DEBUG_CHECK(stride2(lhs) == 1 || size2(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    assign_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    assign_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 assign_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::dgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		  1, data(m2), leading_dimension(m2), 
		  data(m1), leading_dimension(m1),
		  0, data(lhs), stride1(lhs));
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::ContiguousMatrix<double, ColMajor, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, 
                     Multiplication<double, double> f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride1(lhs) == 1 || size1(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    assign_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    assign_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 assign_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::dgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		  1, data(m1), leading_dimension(m1), 
		  data(m2), leading_dimension(m2),
		  0, data(lhs), stride2(lhs));
   }
};

//
// real add
//

template <typename LHS, typename M1, typename M2,
	  typename LHSOrient, typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::StrideMatrix<double, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<double, double> f)
   {
      if (!is_blas_matrix(lhs))
      {
	 Matrix<double, RowMajor> TempL(size1(m1), size2(m2));
	 if (!is_blas_matrix(m1))
	 {
	    Matrix<double, RowMajor> Temp1(m1);
	    if (is_blas_matrix(m2))
	    {
	       assign_product2(TempL, Temp1, m2);
	    }
	    else
	    {
	       Matrix<double, ColMajor> Temp2(m2);
	       assign_product2(TempL, Temp1, Temp2);
	    }
	    add(lhs, TempL);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    add_product2(lhs, Temp1, m2);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    add_product2(lhs, Temp1, Temp2);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 add_product2(lhs, m1, Temp2);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(lhs) && is_blas_matrix(m1) && is_blas_matrix(m2));
      if (blas_trans_col(lhs) == 'T')
      {
	 BLAS::dgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		     1, data(m2), leading_dimension(m2), 
		     data(m1), leading_dimension(m1),
		     1, data(lhs), leading_dimension(lhs));
      }
      else
      {
	 BLAS::dgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		     1, data(m1), leading_dimension(m1), 
		     data(m2), leading_dimension(m2),
		     1, data(lhs), leading_dimension(lhs));
      }
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::ContiguousMatrix<double, RowMajor, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<double, double> f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride2(lhs) == 1 || size2(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    add_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    add_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 add_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::dgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		  1, data(m2), leading_dimension(m2), 
		  data(m1), leading_dimension(m1),
		  1, data(lhs), stride1(lhs));
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::ContiguousMatrix<double, ColMajor, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, 
                     Multiplication<double, double> f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride1(lhs) == 1 || size1(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    add_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    add_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 add_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::dgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		  1, data(m1), leading_dimension(m1), 
		  data(m2), leading_dimension(m2),
		  1, data(lhs), stride2(lhs));
   }
};

//
// real subtract
//

template <typename LHS, typename M1, typename M2,
	  typename LHSOrient, typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::StrideMatrix<double, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<double, double> f)
   {
      if (!is_blas_matrix(lhs))
      {
	 Matrix<double, RowMajor> TempL(size1(m1), size2(m2));
	 if (!is_blas_matrix(m1))
	 {
	    Matrix<double, RowMajor> Temp1(m1);
	    if (is_blas_matrix(m2))
	    {
	       assign_product2(TempL, Temp1, m2);
	    }
	    else
	    {
	       Matrix<double, ColMajor> Temp2(m2);
	       assign_product2(TempL, Temp1, Temp2);
	    }
	    subtract(lhs, TempL);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    subtract_product2(lhs, Temp1, m2);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    subtract_product2(lhs, Temp1, Temp2);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 subtract_product2(lhs, m1, Temp2);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(lhs) && is_blas_matrix(m1) && is_blas_matrix(m2));
      if (blas_trans_col(lhs) == 'T')
      {
	 BLAS::dgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		     -1, data(m2), leading_dimension(m2), 
		     data(m1), leading_dimension(m1),
		     1, data(lhs), leading_dimension(lhs));
      }
      else
      {
	 BLAS::dgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		     -1, data(m1), leading_dimension(m1), 
		     data(m2), leading_dimension(m2),
		     1, data(lhs), leading_dimension(lhs));
      }
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::ContiguousMatrix<double, RowMajor, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<double, double> f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride2(lhs) == 1 || size2(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    subtract_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    subtract_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 subtract_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::dgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		  -1, data(m2), leading_dimension(m2), 
		  data(m1), leading_dimension(m1),
		  1, data(lhs), stride1(lhs));
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Multiplication<double, double>,
		      Concepts::ContiguousMatrix<double, ColMajor, LHSi>,
		      Concepts::StrideMatrix<double, M1Orient, M1i>,
		      Concepts::StrideMatrix<double, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, 
                     Multiplication<double, double> f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride1(lhs) == 1 || size1(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<double, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    subtract_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<double, ColMajor> Temp2(m2);
	    subtract_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<double, ColMajor> Temp2(m2);
	 subtract_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::dgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		  -1, data(m1), leading_dimension(m1), 
		  data(m2), leading_dimension(m2),
		  1, data(lhs), stride2(lhs));
   }
};

//
// =====================================================================================
//

//
// complex assign
//

template <typename LHS, typename M1, typename M2,
	  typename LHSOrient, typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::StrideMatrix<std::complex<double>, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      if (!is_blas_matrix(lhs))
      {
	 Matrix<std::complex<double>, RowMajor> TempL(size1(m1), size2(m2));
	 if (!is_blas_matrix(m1))
	 {
	    Matrix<std::complex<double>, RowMajor> Temp1(m1);
	    if (is_blas_matrix(m2))
	    {
	       assign_product2(TempL, Temp1, m2);
	    }
	    else
	    {
	       Matrix<std::complex<double>, ColMajor> Temp2(m2);
	       assign_product2(TempL, Temp1, Temp2);
	    }
	    assign(lhs, TempL);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    assign_product2(lhs, Temp1, m2);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    assign_product2(lhs, Temp1, Temp2);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 assign_product2(lhs, m1, Temp2);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(lhs) && is_blas_matrix(m1) && is_blas_matrix(m2));
      if (blas_trans_col(lhs) == 'T')
      {
	 BLAS::zgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		     1, data(m2), leading_dimension(m2), 
		     data(m1), leading_dimension(m1),
		     0, data(lhs), leading_dimension(lhs));
      }
      else
      {
	 BLAS::zgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		     1, data(m1), leading_dimension(m1), 
		     data(m2), leading_dimension(m2),
		     0, data(lhs), leading_dimension(lhs));
      }
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::ContiguousMatrix<std::complex<double>, RowMajor, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<std::complex<double>, std::complex<double>  > f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      //      DEBUG_CHECK_EQUAL(stride2(lhs), 1);
      DEBUG_CHECK(stride2(lhs) == 1 || size2(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    assign_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    assign_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 assign_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::zgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		  1, data(m2), leading_dimension(m2), 
		  data(m1), leading_dimension(m1),
		  0, data(lhs), stride1(lhs));
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AssignProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::ContiguousMatrix<std::complex<double>, ColMajor, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, 
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride1(lhs) == 1 || size1(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    assign_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    assign_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 assign_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::zgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		  1, data(m1), leading_dimension(m1), 
		  data(m2), leading_dimension(m2),
		  0, data(lhs), stride2(lhs));
   }
};




//
// complex add
//

template <typename LHS, typename M1, typename M2,
	  typename LHSOrient, typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::StrideMatrix<std::complex<double>, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      if (!is_blas_matrix(lhs))
      {
	 Matrix<std::complex<double>, RowMajor> TempL(size1(m1), size2(m2));
	 if (!is_blas_matrix(m1))
	 {
	    Matrix<std::complex<double>, RowMajor> Temp1(m1);
	    if (is_blas_matrix(m2))
	    {
	       assign_product2(TempL, Temp1, m2);
	    }
	    else
	    {
	       Matrix<std::complex<double>, ColMajor> Temp2(m2);
	       assign_product2(TempL, Temp1, Temp2);
	    }
	    add(lhs, TempL);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    add_product2(lhs, Temp1, m2);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    add_product2(lhs, Temp1, Temp2);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 add_product2(lhs, m1, Temp2);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(lhs) && is_blas_matrix(m1) && is_blas_matrix(m2));
      if (blas_trans_col(lhs) == 'T')
      {
	 BLAS::zgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		     1, data(m2), leading_dimension(m2), 
		     data(m1), leading_dimension(m1),
		     1, data(lhs), leading_dimension(lhs));
      }
      else
      {
	 BLAS::zgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		     1, data(m1), leading_dimension(m1), 
		     data(m2), leading_dimension(m2),
		     1, data(lhs), leading_dimension(lhs));
      }
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::ContiguousMatrix<std::complex<double>, RowMajor, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride2(lhs) == 1 || size2(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    add_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    add_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 add_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::zgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		  1, data(m2), leading_dimension(m2), 
		  data(m1), leading_dimension(m1),
		  1, data(lhs), stride1(lhs));
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct AddProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::ContiguousMatrix<std::complex<double>, ColMajor, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, 
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride1(lhs) == 1 || size1(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    add_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    add_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 add_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::zgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		  1, data(m1), leading_dimension(m1), 
		  data(m2), leading_dimension(m2),
		  1, data(lhs), stride2(lhs));
   }
};

//
// complex subtract
//

template <typename LHS, typename M1, typename M2,
	  typename LHSOrient, typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::StrideMatrix<std::complex<double>, LHSOrient, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      if (!is_blas_matrix(lhs))
      {
	 Matrix<std::complex<double>, RowMajor> TempL(size1(m1), size2(m2));
	 if (!is_blas_matrix(m1))
	 {
	    Matrix<std::complex<double>, RowMajor> Temp1(m1);
	    if (is_blas_matrix(m2))
	    {
	       assign_product2(TempL, Temp1, m2);
	    }
	    else
	    {
	       Matrix<std::complex<double>, ColMajor> Temp2(m2);
	       assign_product2(TempL, Temp1, Temp2);
	    }
	    subtract(lhs, TempL);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    subtract_product2(lhs, Temp1, m2);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    subtract_product2(lhs, Temp1, Temp2);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 subtract_product2(lhs, m1, Temp2);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(lhs) && is_blas_matrix(m1) && is_blas_matrix(m2));
      if (blas_trans_col(lhs) == 'T')
      {
	 BLAS::zgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		     -1, data(m2), leading_dimension(m2), 
		     data(m1), leading_dimension(m1),
		     1, data(lhs), leading_dimension(lhs));
      }
      else
      {
	 BLAS::zgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		     -1, data(m1), leading_dimension(m1), 
		     data(m2), leading_dimension(m2),
		     1, data(lhs), leading_dimension(lhs));
      }
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::ContiguousMatrix<std::complex<double>, RowMajor, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2,
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride2(lhs) == 1 || size2(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    subtract_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    subtract_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 subtract_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::zgemm(blas_trans_row(m2), blas_trans_row(m1), size2(m2), size1(m1), size2(m1),
		  -1, data(m2), leading_dimension(m2), 
		  data(m1), leading_dimension(m1),
		  1, data(lhs), stride1(lhs));
   }
};

template <typename LHS, typename M1, typename M2,
	  typename LHSi,
	  typename M1Orient, typename M1i,
	  typename M2Orient, typename M2i>
struct SubtractProduct2<LHS, M1, M2, Multiplication<std::complex<double>, std::complex<double> >,
		      Concepts::ContiguousMatrix<std::complex<double>, ColMajor, LHSi>,
		      Concepts::StrideMatrix<std::complex<double>, M1Orient, M1i>,
		      Concepts::StrideMatrix<std::complex<double>, M2Orient, M2i>>
{
   static void apply(LHS& lhs, M1 const& m1, M2 const& m2, 
                     Multiplication<std::complex<double>, std::complex<double> > f)
   {
      DEBUG_CHECK(is_blas_matrix(lhs));
      DEBUG_CHECK(stride1(lhs) == 1 || size1(lhs) == 1);
      if (!is_blas_matrix(m1))
      {
	 Matrix<std::complex<double>, RowMajor> Temp1(m1);
	 if (is_blas_matrix(m2))
	 {
	    subtract_product2(lhs, Temp1, m2, f);
	 }
	 else
	 {
	    Matrix<std::complex<double>, ColMajor> Temp2(m2);
	    subtract_product2(lhs, Temp1, Temp2, f);
	 }
	 return;
      }
      // else
      if (!is_blas_matrix(m2))
      {
	 Matrix<std::complex<double>, ColMajor> Temp2(m2);
	 subtract_product2(lhs, m1, Temp2, f);
         return;
      }
      // else
      DEBUG_CHECK(is_blas_matrix(m1) && is_blas_matrix(m2));
      BLAS::zgemm(blas_trans_col(m1), blas_trans_col(m2), size1(m1), size2(m2), size2(m1),
		  -1, data(m1), leading_dimension(m1), 
		  data(m2), leading_dimension(m2),
		  1, data(lhs), stride2(lhs));
   }
};












} // namespace LinearAlgebra

#endif
