// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/matrixproductold.h
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



//
// assign_product3
//

template <typename Scalar, typename Derived,
          typename S1, typename D1, typename S2, typename D2, typename S3, typename D3>
void
assign_product3(MatrixExpression<Scalar, Derived>& lhs,
                MatrixConstExpression<S1, D1> const& r1,
                MatrixConstExpression<S2, D2> const& r2,
                MatrixConstExpression<S3, D3> const& r3)
{
   // choose an order of evaluation that minimizes the number of operations
   if (r1.size1() * r2.size2() * (r2.size1() + r3.size2())
       < (r1.size1() + r2.size2()) * r2.size1() * r3.size2())
   {
      typedef typename BinaryOperator<Multiplication, S1, S2>::value_type S1S2_value_type;
      TempMatrix<S1S2_value_type> Temp12(r1.size1(), r2.size2());
      assign_product2(Temp12, r1.as_derived(), r2.as_derived());
      assign_product2(lhs.as_derived(), Temp12, r3.as_derived());
   }

   {
      typedef typename BinaryOperator<Multiplication, S2, S3>::value_type S2S3_value_type;
      TempMatrix<S2S3_value_type> Temp23(r2.size1(), r3.size2());
      assign_product2(Temp23, r2.as_derived(), r3.as_derived());
      assign_product2(lhs.as_derived(), r1.as_derived(), Temp23);
   }
}

//
// assign_product2
//

#if !defined(LINEARALGEBRA_NO_TEMP_SPECIALIZATION)

template <typename Scalar, typename Derived, typename S1, typename D1, typename S2, typename D2>
void
assign_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixConstExpression<S1, D1> const& r1,
                MatrixConstExpression<S2, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   typedef typename MatrixExpression<Scalar, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                          iiterator;
   typedef typename D1::const_iterator1                          M1iterator1;
   typedef typename D2::const_iterator2                          M2iterator2;

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

#else

template <typename Scalar, typename Derived, typename S1, typename D1, typename S2, typename D2>
void
assign_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixConstSlice<S1, D1> const& r1,
                MatrixConstSlice<S2, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   typedef typename MatrixExpression<Scalar, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                          iiterator;
   typedef typename D1::const_iterator1                          M1iterator1;
   typedef typename D2::const_iterator2                          M2iterator2;

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

#endif

template <typename Scalar, typename Derived, typename S1, typename D1, typename S2, typename D2>
void
assign_product2(MatrixExpression<Scalar, Derived>& lhs,
                GenericMatrix<S1, D1> const& r1,
                GenericMatrix<S2, D2> const& r2)
{
   TempMatrix<S1> Temp1(r1);
   TempMatrix<S2> Temp2(r2);
   assign_product2(lhs, Temp1, Temp2);
}

template <typename Scalar, typename Derived, typename M1, typename M2, typename S, typename D>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixProductExpression<M1, M2> const& r1,
                MatrixConstExpression<S, D> const& r2)
{
   assign_product3(lhs.as_derived(), r1.matrix1().eval_expr(),
                   r1.matrix2().eval_expr(), r2.as_derived());
}

template <typename Scalar, typename Derived, typename M1, typename M2, typename S, typename D>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixConstExpression<S, D> const& r1,
                MatrixProductExpression<M1, M2> const& r2)
{
   assign_product3(lhs.as_derived(), r1.as_derived(),
                   r2.matrix1().eval_expr(), r2.matrix2().eval_expr());
}

template <typename Scalar, typename Derived, typename M1, typename M2, typename M3, typename M4>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixProductExpression<M1, M2> const& r1,
                MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   assign_product2(lhs.as_derived(), Temp1, Temp2);
}

//
// assign_scaled_product3
//

template <typename Scalar, typename Derived, typename S,
          typename S1, typename D1, typename S2, typename D2, typename S3, typename D3>
void
assign_scaled_product3(MatrixExpression<Scalar, Derived>& lhs,
                       S const& x,
                       GenericMatrix<S1, D1> const& r1,
                       GenericMatrix<S2, D2> const& r2,
                       GenericMatrix<S3, D3> const& r3)
{
   // choose an order of evaluation that minimizes the number of operations
   if (r1.size1() * r2.size2() * (r2.size1() + r3.size2())
       < (r1.size1() + r2.size2()) * r2.size1() * r3.size2())
   {
      typedef typename BinaryOperator<Multiplication, S1, S2>::value_type S1S2_value_type;
      TempMatrix<S1S2_value_type> Temp12(r1.size1(), r2.size2());
      assign_product2(Temp12, r1.as_derived(), r2.as_derived());
      assign_scaled_product2(lhs.as_derived(), x, Temp12, r3.as_derived());
   }

   {
      typedef typename BinaryOperator<Multiplication, S2, S3>::value_type S2S3_value_type;
      TempMatrix<S2S3_value_type> Temp23(r2.size1(), r3.size2());
      assign_product2(Temp23, r2.as_derived(), r3.as_derived());
      assign_scaled_product2(lhs.as_derived(), x, r1.as_derived(), Temp23);
   }
}

//
// assign_scaled_product2
//

template <typename Scalar, typename Derived,
          typename S, typename S1, typename D1, typename S2, typename D2>
void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                       S const& x,
                       MatrixConstExpression<S1, D1> const& r1,
                       MatrixConstExpression<S2, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   typedef typename MatrixExpression<Scalar, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                          iiterator;
   typedef typename D1::const_iterator1                          M1iterator1;
   typedef typename D2::const_iterator2                          M2iterator2;

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

template <typename Scalar, typename Derived, typename Sc,
          typename M1, typename M2, typename S, typename D>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                       Sc const& x,
                       MatrixProductExpression<M1, M2> const& r1,
                       MatrixConstExpression<S, D> const& r2)
{
   assign_scaled_product3(lhs.as_derived(), x, r1.matrix1(), r1.matrix2(), r2.as_derived());
}

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename S, typename D>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                       Sc const& x,
                       MatrixConstExpression<S, D> const& r1,
                       MatrixProductExpression<M1, M2> const& r2)
{
   assign_scaled_product3(lhs.as_derived(), x, r1.as_derived(), r2.matrix1(), r2.matrix2());
}

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename M3, typename M4>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                       Sc const& x,
                       MatrixProductExpression<M1, M2> const& r1,
                       MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   assign_scaled_product2(lhs.as_derived(), x, Temp1, Temp2);
}

//
// add_product3
//

template <typename Scalar, typename Derived,
          typename S1, typename D1, typename S2, typename D2, typename S3, typename D3>
void
add_product3(MatrixExpression<Scalar, Derived>& lhs,
             GenericMatrix<S1, D1> const& r1,
             GenericMatrix<S2, D2> const& r2,
             GenericMatrix<S3, D3> const& r3)
{
   // choose an order of evaluation that minimizes the number of operations
   if (r1.size1() * r2.size2() * (r2.size1() + r3.size2())
       < (r1.size1() + r2.size2()) * r2.size1() * r3.size2())
   {
      typedef typename BinaryOperator<Multiplication, S1, S2>::value_type S1S2_value_type;
      TempMatrix<S1S2_value_type> Temp12(r1.size1(), r2.size2());
      assign_product2(Temp12, r1.as_derived(), r2.as_derived());
      add_product2(lhs.as_derived(), Temp12, r3.as_derived());
   }

   {
      typedef typename BinaryOperator<Multiplication, S2, S3>::value_type S2S3_value_type;
      TempMatrix<S2S3_value_type> Temp23(r2.size1(), r3.size2());
      assign_product2(Temp23, r2.as_derived(), r3.as_derived());
      add_product2(lhs.as_derived(), r1.as_derived(), Temp23);
   }
}

//
// add_product2
//

template <typename Scalar, typename Derived, typename S1, typename D1, typename S2, typename D2>
void
add_product2(MatrixExpression<Scalar, Derived>& lhs,
             MatrixConstExpression<S1, D1> const& r1,
             MatrixConstExpression<S2, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   typedef typename MatrixExpression<Scalar, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                          iiterator;
   typedef typename D1::const_iterator1                          M1iterator1;
   typedef typename D2::const_iterator2                          M2iterator2;

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

template <typename Scalar, typename Derived, typename M1, typename M2, typename S, typename D>
inline void
add_product2(MatrixExpression<Scalar, Derived>& lhs,
             MatrixProductExpression<M1, M2> const& r1,
             MatrixConstExpression<S, D> const& r2)
{
   add_product3(lhs.as_derived(), r1.matrix1(), r1.matrix2(), r2.as_derived());
}

template <typename Scalar, typename Derived, typename M1, typename M2, typename S, typename D>
inline void
add_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixConstExpression<S, D> const& r1,
                MatrixProductExpression<M1, M2> const& r2)
{
   add_product3(lhs.as_derived(), r1.as_derived(), r2.matrix1(), r2.matrix2());
}

template <typename Scalar, typename Derived, typename M1, typename M2, typename M3, typename M4>
inline void
add_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixProductExpression<M1, M2> const& r1,
                MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   add_product2(lhs.as_derived(), Temp1, Temp2);
}

//
// add_scaled_product3
//

template <typename Scalar, typename Derived, typename S,
          typename S1, typename D1, typename S2, typename D2, typename S3, typename D3>
void
add_scaled_product3(MatrixExpression<Scalar, Derived>& lhs,
                    S const& x,
                    GenericMatrix<S1, D1> const& r1,
                    GenericMatrix<S2, D2> const& r2,
                    GenericMatrix<S3, D3> const& r3)
{
   // choose an order of evaluation that minimizes the number of operations
   if (r1.size1() * r2.size2() * (r2.size1() + r3.size2())
       < (r1.size1() + r2.size2()) * r2.size1() * r3.size2())
   {
      typedef typename BinaryOperator<Multiplication, S1, S2>::value_type S1S2_value_type;
      TempMatrix<S1S2_value_type> Temp12(r1.size1(), r2.size2());
      assign_product2(Temp12, r1.as_derived(), r2.as_derived());
      add_scaled_product2(lhs.as_derived(), x, Temp12, r3.as_derived());
   }

   {
      typedef typename BinaryOperator<Multiplication, S2, S3>::value_type S2S3_value_type;
      TempMatrix<S2S3_value_type> Temp23(r2.size1(), r3.size2());
      assign_product2(Temp23, r2.as_derived(), r3.as_derived());
      add_scaled_product2(lhs.as_derived(), x, r1.as_derived(), Temp23);
   }
}

//
// add_scaled_product2
//

template <typename Scalar, typename Derived,
          typename S, typename S1, typename D1, typename S2, typename D2>
void
add_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    S const& x,
                    MatrixConstExpression<S1, D1> const& r1,
                    MatrixConstExpression<S2, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   typedef typename MatrixExpression<Scalar, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                          iiterator;
   typedef typename D1::const_iterator1                          M1iterator1;
   typedef typename D2::const_iterator2                          M2iterator2;

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

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename S, typename D>
inline void
add_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    Sc const& x,
                    MatrixProductExpression<M1, M2> const& r1,
                    MatrixConstExpression<S, D> const& r2)
{
   add_scaled_product3(lhs.as_derived(), x, r1.matrix1(), r1.matrix2(), r2.as_derived());
}

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename S, typename D>
inline void
add_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    Sc const& x,
                    MatrixConstExpression<S, D> const& r1,
                    MatrixProductExpression<M1, M2> const& r2)
{
   add_scaled_product3(lhs.as_derived(), x, r1.as_derived(), r2.matrix1(), r2.matrix2());
}

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename M3, typename M4>
inline void
add_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    Sc const& x,
                    MatrixProductExpression<M1, M2> const& r1,
                    MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   add_scaled_product2(lhs.as_derived(), x, Temp1, Temp2);
}

//
// sub_product3
//

template <typename Scalar, typename Derived,
          typename S1, typename D1, typename S2, typename D2, typename S3, typename D3>
void
sub_product3(MatrixExpression<Scalar, Derived>& lhs,
             GenericMatrix<S1, D1> const& r1,
             GenericMatrix<S2, D2> const& r2,
             GenericMatrix<S3, D3> const& r3)
{
   // choose an order of evaluation that minimizes the number of operations
   if (r1.size1() * r2.size2() * (r2.size1() + r3.size2())
       < (r1.size1() + r2.size2()) * r2.size1() * r3.size2())
   {
      typedef typename BinaryOperator<Multiplication, S1, S2>::value_type S1S2_value_type;
      TempMatrix<S1S2_value_type> Temp12(r1.size1(), r2.size2());
      assign_product2(Temp12, r1.as_derived(), r2.as_derived());
      sub_product2(lhs.as_derived(), Temp12, r3.as_derived());
   }

   {
      typedef typename BinaryOperator<Multiplication, S2, S3>::value_type S2S3_value_type;
      TempMatrix<S2S3_value_type> Temp23(r2.size1(), r3.size2());
      assign_product2(Temp23, r2.as_derived(), r3.as_derived());
      sub_product2(lhs.as_derived(), r1.as_derived(), Temp23);
   }
}

//
// sub_product2
//

template <typename Scalar, typename Derived, typename S1, typename D1, typename S2, typename D2>
void
sub_product2(MatrixExpression<Scalar, Derived>& lhs,
             MatrixConstExpression<S1, D1> const& r1,
             MatrixConstExpression<S2, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   typedef typename MatrixExpression<Scalar, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                          iiterator;
   typedef typename D1::const_iterator1                          M1iterator1;
   typedef typename D2::const_iterator2                          M2iterator2;

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

template <typename Scalar, typename Derived, typename M1, typename M2, typename S, typename D>
inline void
sub_product2(MatrixExpression<Scalar, Derived>& lhs,
             MatrixProductExpression<M1, M2> const& r1,
             MatrixConstExpression<S, D> const& r2)
{
   sub_product3(lhs.as_derived(), r1.matrix1(), r1.matrix2(), r2.as_derived());
}

template <typename Scalar, typename Derived, typename M1, typename M2, typename S, typename D>
inline void
sub_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixConstExpression<S, D> const& r1,
                MatrixProductExpression<M1, M2> const& r2)
{
   sub_product3(lhs.as_derived(), r1.as_derived(), r2.matrix1(), r2.matrix2());
}

template <typename Scalar, typename Derived, typename M1, typename M2, typename M3, typename M4>
inline void
sub_product2(MatrixExpression<Scalar, Derived>& lhs,
                MatrixProductExpression<M1, M2> const& r1,
                MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   sub_product2(lhs.as_derived(), Temp1, Temp2);
}

//
// sub_scaled_product3
//

template <typename Scalar, typename Derived, typename S,
          typename S1, typename D1, typename S2, typename D2, typename S3, typename D3>
void
sub_scaled_product3(MatrixExpression<Scalar, Derived>& lhs,
                    S const& x,
                    GenericMatrix<S1, D1> const& r1,
                    GenericMatrix<S2, D2> const& r2,
                    GenericMatrix<S3, D3> const& r3)
{
   // choose an order of evaluation that minimizes the number of operations
   if (r1.size1() * r2.size2() * (r2.size1() + r3.size2())
       < (r1.size1() + r2.size2()) * r2.size1() * r3.size2())
   {
      typedef typename BinaryOperator<Multiplication, S1, S2>::value_type S1S2_value_type;
      TempMatrix<S1S2_value_type> Temp12(r1.size1(), r2.size2());
      assign_product2(Temp12, r1.as_derived(), r2.as_derived());
      sub_scaled_product2(lhs.as_derived(), x, Temp12, r3.as_derived());
   }

   {
      typedef typename BinaryOperator<Multiplication, S2, S3>::value_type S2S3_value_type;
      TempMatrix<S2S3_value_type> Temp23(r2.size1(), r3.size2());
      assign_product2(Temp23, r2.as_derived(), r3.as_derived());
      sub_scaled_product2(lhs.as_derived(), x, r1.as_derived(), Temp23);
   }
}

//
// sub_scaled_product2
//

template <typename Scalar, typename Derived,
          typename S, typename S1, typename D1, typename S2, typename D2>
void
sub_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    S const& x,
                    MatrixConstExpression<S1, D1> const& r1,
                    MatrixConstExpression<S2, D2> const& r2)
{
   PRECONDITION(lhs.size1() == r1.size1())(lhs.size1())(r1.size1());
   PRECONDITION(lhs.size2() == r2.size2())(lhs.size2())(r2.size2());
   PRECONDITION(r1.size2() == r2.size1())(r1.size2())(r2.size1());

   typedef typename MatrixExpression<Scalar, Derived>::iterator1 oiterator;
   typedef typename oiterator::iterator                          iiterator;
   typedef typename D1::const_iterator1                          M1iterator1;
   typedef typename D2::const_iterator2                          M2iterator2;

   oiterator End = lhs.end1();
   M1iterator1 M1I = r1.begin1();
   for (oiterator I = lhs.begin1(); I != End; ++I, ++M1I)
   {
      iiterator JEnd = I.end();
      M2iterator2 M2J = r2.begin2();
      for (iiterator J = I.begin(); J != JEnd; ++J, ++M2J)
      {
         *J -= x * inner_prod(*M1I, *M2J);
      }
   }
}

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename S, typename D>
inline void
sub_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    Sc const& x,
                    MatrixProductExpression<M1, M2> const& r1,
                    MatrixConstExpression<S, D> const& r2)
{
   sub_scaled_product3(lhs.as_derived(), x, r1.matrix1(), r1.matrix2(), r2.as_derived());
}

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename S, typename D>
inline void
sub_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    Sc const& x,
                    MatrixConstExpression<S, D> const& r1,
                    MatrixProductExpression<M1, M2> const& r2)
{
   sub_scaled_product3(lhs.as_derived(), x, r1.as_derived(), r2.matrix1(), r2.matrix2());
}

template <typename Scalar, typename Derived,
          typename Sc, typename M1, typename M2, typename M3, typename M4>
inline void
sub_scaled_product2(MatrixExpression<Scalar, Derived>& lhs,
                    Sc const& x,
                    MatrixProductExpression<M1, M2> const& r1,
                    MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   sub_scaled_product2(lhs.as_derived(), x, Temp1, Temp2);
}
