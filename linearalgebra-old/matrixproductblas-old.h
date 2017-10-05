// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrixproductblas-old.h
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


#if !defined(MATRIXPRODUCTBLAS_H_DSHJFW5Y789FY78943FOH)
#define MATRIXPRODUCTBLAS_H_DSHJFW5Y789FY78943FOH

namespace LinearAlgebra
{

//
// some utility functions to determine the TRANS parameter for BLAS
//

template <typename Scalar, typename Derived>
inline
bool is_blas(MatrixConstSlice<Scalar, Derived> const& M)
{
   return M.stride1() == 1 || M.stride2() == 1;
}

template <typename Scalar, typename Derived>
inline
char blas_trans_col(MatrixConstSlice<Scalar, Derived> const& M)
{
   // note: it is important that we are consistent with comparing stride2() first
   // in blas_trans_col, blas_trans_row and leading_dimension.  Otherwise
   // there is a problem with Mx1 and 1xM matrices.
   if (M.stride2() == 1) return 'T';
   if (M.stride1() == 1) return 'N';
   return 'X';
}

template <typename Scalar, typename Derived>
inline
char blas_trans_row(MatrixConstSlice<Scalar, Derived> const& M)
{
   if (M.stride2() == 1) return 'N';
   if (M.stride1() == 1) return 'T';
   return 'X';
}

template <typename Scalar, typename Derived>
inline
size_t leading_dimension(MatrixConstSlice<Scalar, Derived> const& M)
{
   if (M.stride2() == 1) return M.stride1();
   return M.stride2();
}

//
// real
//

template <typename Derived, typename D1, typename D2>
void
assign_product2(MatrixSlice<double, Derived>& lhs,
                MatrixConstSlice<double, D1> const& r1,
                MatrixConstSlice<double, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   if (is_blas(lhs) && is_blas(r1) && is_blas(r2))
   {
      if (blas_trans_col(lhs) == 'T')
      {
         BLAS::dgemm(blas_trans_row(r2), blas_trans_row(r1), r2.size2(), r1.size1(), r1.size2(),
                     1, r2.data(), leading_dimension(r2),
                     r1.data(), leading_dimension(r1),
                     0, lhs.data(), leading_dimension(lhs));
      }
      else
      {
         BLAS::dgemm(blas_trans_col(r1), blas_trans_col(r2), r1.size1(), r2.size2(), r1.size2(),
                     1, r1.data(), leading_dimension(r1),
                     r2.data(), leading_dimension(r2),
                     0, lhs.data(), leading_dimension(lhs));
      }
      return;
   }

   // else do it the slow way

   typedef typename MatrixSlice<double, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                     iiterator;
   typedef typename D1::const_iterator1                     M1iterator1;
   typedef typename D2::const_iterator2                     M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J = inner_prod(*M1I, *M2J);
      }
   }
}

template <typename Derived, typename D1, typename D2>
void
add_product2(MatrixSlice<double, Derived>& lhs,
             MatrixConstSlice<double, D1> const& r1,
             MatrixConstSlice<double, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   if (is_blas(lhs) && is_blas(r1) && is_blas(r2))
   {
      if (blas_trans_col(lhs) == 'T')
      {
         BLAS::dgemm(blas_trans_row(r2), blas_trans_row(r1), r2.size2(), r1.size1(), r1.size2(),
                     1, r2.data(), leading_dimension(r2),
                     r1.data(), leading_dimension(r1),
                     1, lhs.data(), leading_dimension(lhs));
      }
      else
      {
         BLAS::dgemm(blas_trans_col(r1), blas_trans_col(r2), r1.size1(), r2.size2(), r1.size2(),
                     1, r1.data(), leading_dimension(r1),
                     r2.data(), leading_dimension(r2),
                     1, lhs.data(), leading_dimension(lhs));
      }
      return;
   }

   // else do it the slow way

   typedef typename MatrixSlice<double, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                     iiterator;
   typedef typename D1::const_iterator1                     M1iterator1;
   typedef typename D2::const_iterator2                     M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J += inner_prod(*M1I, *M2J);
      }
   }
}

template <typename Derived, typename D1, typename D2>
void
sub_product2(MatrixSlice<double, Derived>& lhs,
             MatrixConstSlice<double, D1> const& r1,
             MatrixConstSlice<double, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   if (is_blas(lhs) && is_blas(r1) && is_blas(r2))
   {
      if (blas_trans_col(lhs) == 'T')
      {
         BLAS::dgemm(blas_trans_row(r2), blas_trans_row(r1), r2.size2(), r1.size1(), r1.size2(),
                     -1, r2.data(), leading_dimension(r2),
                     r1.data(), leading_dimension(r1),
                     1, lhs.data(), leading_dimension(lhs));
      }
      else
      {
         BLAS::dgemm(blas_trans_col(r1), blas_trans_col(r2), r1.size1(), r2.size2(), r1.size2(),
                     -1, r1.data(), leading_dimension(r1),
                     r2.data(), leading_dimension(r2),
                     1, lhs.data(), leading_dimension(lhs));
      }
      return;
   }

   // else do it the slow way

   typedef typename MatrixSlice<double, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                     iiterator;
   typedef typename D1::const_iterator1                     M1iterator1;
   typedef typename D2::const_iterator2                     M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J -= inner_prod(*M1I, *M2J);
      }
   }
}

template <typename Derived, typename Scalar, typename D1, typename D2>
void
assign_scaled_product2(MatrixSlice<double, Derived>& lhs,
                       Scalar x,
                       MatrixConstSlice<double, D1> const& r1,
                       MatrixConstSlice<double, D2> const& r2,
                       typename boost::enable_if<
                       boost::is_convertible<Scalar, double> >::type* dummy = 0)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   if (is_blas(lhs) && is_blas(r1) && is_blas(r2))
   {
      if (blas_trans_col(lhs) == 'T')
      {
         BLAS::dgemm(blas_trans_row(r2), blas_trans_row(r1), r2.size2(), r1.size1(), r1.size2(),
                     x, r2.data(), leading_dimension(r2),
                     r1.data(), leading_dimension(r1),
                     0, lhs.data(), leading_dimension(lhs));
      }
      else
      {
         BLAS::dgemm(blas_trans_col(r1), blas_trans_col(r2), r1.size1(), r2.size2(), r1.size2(),
                     x, r1.data(), leading_dimension(r1),
                     r2.data(), leading_dimension(r2),
                     0, lhs.data(), leading_dimension(lhs));
      }
      return;
   }

   // else do it the slow way

   typedef typename MatrixSlice<double, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                     iiterator;
   typedef typename D1::const_iterator1                     M1iterator1;
   typedef typename D2::const_iterator2                     M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J = x * inner_prod(*M1I, *M2J);
      }
   }
}


template <typename Derived, typename Scalar, typename D1, typename D2>
void
add_scaled_product2(MatrixSlice<double, Derived>& lhs,
                    Scalar x,
                    MatrixConstSlice<double, D1> const& r1,
                    MatrixConstSlice<double, D2> const& r2,
                    typename boost::enable_if<boost::is_convertible<Scalar, double> >::type* dummy = 0)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   if (is_blas(lhs) && is_blas(r1) && is_blas(r2))
   {
      if (blas_trans_col(lhs) == 'T')
      {
         BLAS::dgemm(blas_trans_row(r2), blas_trans_row(r1), r2.size2(), r1.size1(), r1.size2(),
                     x, r2.data(), leading_dimension(r2),
                     r1.data(), leading_dimension(r1),
                     1, lhs.data(), leading_dimension(lhs));
      }
      else
      {
         BLAS::dgemm(blas_trans_col(r1), blas_trans_col(r2), r1.size1(), r2.size2(), r1.size2(),
                     x, r1.data(), leading_dimension(r1),
                     r2.data(), leading_dimension(r2),
                     1, lhs.data(), leading_dimension(lhs));
      }
      return;
   }

   // else do it the slow way

   typedef typename MatrixSlice<double, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                     iiterator;
   typedef typename D1::const_iterator1                     M1iterator1;
   typedef typename D2::const_iterator2                     M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J += x * inner_prod(*M1I, *M2J);
      }
   }
}

template <typename Derived, typename Scalar, typename D1, typename D2>
void
sub_scaled_product2(MatrixSlice<double, Derived>& lhs,
                    Scalar x,
                    MatrixConstSlice<double, D1> const& r1,
                    MatrixConstSlice<double, D2> const& r2,
                    typename boost::enable_if<
                      boost::is_convertible<Scalar, double> >::type* dummy = 0)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   if (is_blas(lhs) && is_blas(r1) && is_blas(r2))
   {
      if (blas_trans_col(lhs) == 'T')
      {
         BLAS::dgemm(blas_trans_row(r2), blas_trans_row(r1), r2.size2(), r1.size1(), r1.size2(),
                     -x, r2.data(), leading_dimension(r2),
                     r1.data(), leading_dimension(r1),
                     1, lhs.data(), leading_dimension(lhs));
      }
      else
      {
         BLAS::dgemm(blas_trans_col(r1), blas_trans_col(r2), r1.size1(), r2.size2(), r1.size2(),
                     -x, r1.data(), leading_dimension(r1),
                     r2.data(), leading_dimension(r2),
                     1, lhs.data(), leading_dimension(lhs));
      }
      return;
   }

   // else do it the slow way

   typedef typename MatrixSlice<double, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                     iiterator;
   typedef typename D1::const_iterator1                     M1iterator1;
   typedef typename D2::const_iterator2                     M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J += x * inner_prod(*M1I, *M2J);
      }
   }
}


//
// complex
//

template <typename Derived, typename D1, typename D2>
void
assign_product2(MatrixSlice<std::complex<double>, Derived>& lhs,
                MatrixConstSlice<std::complex<double>, D1> const& r1,
                MatrixConstSlice<std::complex<double>, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   if (is_blas(lhs) && is_blas(r1) && is_blas(r2))
   {
      if (blas_trans_col(lhs) == 'T')
      {
         BLAS::zgemm(blas_trans_row(r2), blas_trans_row(r1), r2.size2(), r1.size1(), r1.size2(),
                     1, r2.data(), leading_dimension(r2),
                     r1.data(), leading_dimension(r1),
                     0, lhs.data(), leading_dimension(lhs));
      }
      else
      {
         BLAS::zgemm(blas_trans_col(r1), blas_trans_col(r2), r1.size1(), r2.size2(), r1.size2(),
                     1, r1.data(), leading_dimension(r1),
                     r2.data(), leading_dimension(r2),
                     0, lhs.data(), leading_dimension(lhs));
      }
      return;
   }

   // else do it the slow way

   typedef typename MatrixSlice<std::complex<double>, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                     iiterator;
   typedef typename D1::const_iterator1                     M1iterator1;
   typedef typename D2::const_iterator2                     M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J = inner_prod(*M1I, *M2J);
      }
   }
}

template <typename Derived, typename D1, typename D2>
void
assign_product2(MatrixSlice<std::complex<double>, Derived>& lhs,
                MatrixConstExpression<std::complex<double>, D1> const& r1,
                MatrixConstExpression<std::complex<double>, D2> const& r2)
{
   TempMatrix<std::complex<double> > x1(r1.as_derived());
   TempMatrix<std::complex<double> > x2(r2.as_derived());
   assign_product2(lhs.as_derived(), x1, x2);
}

} // namespace LinearAlgebra

#endif
