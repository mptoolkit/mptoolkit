// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/matrixproduct-removed.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
#if 0
// FIXME: I don't think these are ever called - prod() unwraps the multiply by scalar itself.
// The only time we don't want to do that is when the multiply doesn't commute/associate properly.
// in which case these overloads are incorrect anyway

// adaptors for ScaledMatrixProductExpression

template <typename Scalar, typename Derived, typename Sc, typename M1, typename M2, typename S, typename D>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs, 
		ScaledMatrixProductExpression<Sc, M1, M2> const& r1,
		MatrixConstExpression<S, D> const& r2)
{
   assign_scaled_product3(lhs.as_derived(), r1.value(), r1.matrix1(), r1.matrix2(), r2.as_derived());
}

template <typename Scalar, typename Derived, typename Sc, typename M1, typename M2, typename S, typename D>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs, 
		MatrixConstExpression<S, D> const& r1,
		ScaledMatrixProductExpression<Sc, M1, M2> const& r2)
{
   assign_product3(lhs.as_derived(), r1.value(), r1.as_derived(), r2.matrix1(), r2.matrix2());
}

template <typename Scalar, typename Derived, typename Sc, typename M1, typename M2, typename M3, typename M4>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs, 
		ScaledMatrixProductExpression<Sc, M1, M2> const& r1,
		MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename ScaledMatrixProductExpression<Sc, M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_scaled_product2(Temp1, r1.value(), r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   assign_product2(lhs.as_derived(), Temp1, Temp2);
}

template <typename Scalar, typename Derived, typename Sc, typename M1, typename M2, typename M3, typename M4>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs, 
		MatrixProductExpression<M1, M2> const& r1,
		ScaledMatrixProductExpression<Sc, M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename ScaledMatrixProductExpression<Sc, M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_scaled_product2(Temp2, r2.value(), r2.matrix1(), r2.matrix2());

   assign_product2(lhs.as_derived(), Temp1, Temp2);
}

template <typename Scalar, typename Derived, typename S1, typename S2,
	  typename M1, typename M2, typename M3, typename M4>
inline void
assign_product2(MatrixExpression<Scalar, Derived>& lhs, 
		ScaledMatrixProductExpression<S1, M1, M2> const& r1,
		ScaledMatrixProductExpression<S2, M3, M4> const& r2)
{
   TempMatrix<typename ScaledMatrixProductExpression<S1, M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_scaled_product2(Temp1, r1.value(), r1.matrix1(), r1.matrix2());

   TempMatrix<typename ScaledMatrixProductExpression<S2, M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_scaled_product2(Temp2, r2.value(), r2.matrix1(), r2.matrix2());

   assign_product2(lhs.as_derived(), Temp1, Temp2);
}

#endif









// adaptors for ScaledMatrixProductExpression

#if 0
// FIXME: again, these overloads should never be used

template <typename Scalar, typename Derived, typename SS, typename Sc, typename M1, typename M2, typename S, typename D>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs, 
		       SS const& x,
		       ScaledMatrixProductExpression<Sc, M1, M2> const& r1,
		       MatrixConstExpression<S, D> const& r2)
{
   assign_scaled_product3(lhs.as_derived(), x * r1.value(), r1.matrix1(), r1.matrix2(), r2.as_derived());
}

template <typename Scalar, typename Derived, typename SS, typename Sc, typename M1, typename M2, typename S, typename D>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs, 
		       SS const& x,
		       MatrixConstExpression<S, D> const& r1,
		       ScaledMatrixProductExpression<Sc, M1, M2> const& r2)
{
   assign_product3(lhs.as_derived(), x* r1.value(), r1.as_derived(), r2.matrix1(), r2.matrix2());
}

template <typename Scalar, typename Derived, typename SS, typename Sc, typename M1, typename M2, typename M3, typename M4>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs, 
		       SS const& x,
		       ScaledMatrixProductExpression<Sc, M1, M2> const& r1,
		       MatrixProductExpression<M3, M4> const& r2)
{
   TempMatrix<typename ScaledMatrixProductExpression<Sc, M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_scaled_product2(Temp1, r1.value(), r1.matrix1(), r1.matrix2());

   TempMatrix<typename MatrixProductExpression<M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_product2(Temp2, r2.matrix1(), r2.matrix2());

   assign_scaled_product2(lhs.as_derived(), x, Temp1, Temp2);
}

template <typename Scalar, typename Derived, typename SS, typename Sc, typename M1, typename M2, typename M3, typename M4>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs, 
		       SS const& x,
		       MatrixProductExpression<M1, M2> const& r1,
		       ScaledMatrixProductExpression<Sc, M3, M4> const& r2)
{
   TempMatrix<typename MatrixProductExpression<M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_product2(Temp1, r1.matrix1(), r1.matrix2());

   TempMatrix<typename ScaledMatrixProductExpression<Sc, M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_scaled_product2(Temp2, r2.value(), r2.matrix1(), r2.matrix2());

   assign_scaled_product2(lhs.as_derived(), x, Temp1, Temp2);
}

template <typename Scalar, typename Derived, typename SS, typename S1, typename S2,
	  typename M1, typename M2, typename M3, typename M4>
inline void
assign_scaled_product2(MatrixExpression<Scalar, Derived>& lhs, 
		       SS const& x,
		       ScaledMatrixProductExpression<S1, M1, M2> const& r1,
		       ScaledMatrixProductExpression<S2, M3, M4> const& r2)
{
   TempMatrix<typename ScaledMatrixProductExpression<S1, M1, M2>::value_type> Temp1(r1.size1(), r1.size2());
   assign_scaled_product2(Temp1, r1.value(), r1.matrix1(), r1.matrix2());

   TempMatrix<typename ScaledMatrixProductExpression<S2, M3, M4>::value_type> Temp2(r2.size1(), r2.size2());
   assign_scaled_product2(Temp2, r2.value(), r2.matrix1(), r2.matrix2());

   assign_scaled_product2(lhs.as_derived(), x, Temp1, Temp2);
}

#endif
