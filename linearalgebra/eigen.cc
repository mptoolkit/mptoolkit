// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/eigen.cc
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

#include "matrix.h"
#include "scalarmatrix.h"
#include "matrix_utility.h"

// randomize signs of vectors in debug mode
#if !defined(NDEBUG) && !defined(RANDOMIZE_VECTORS)
#define RANDOMIZE_VECTORS
#endif

namespace LinearAlgebra
{

namespace Private
{

void LinearSolveSPD(int Size, int Nrhs, double* A, int ldA, double* B, int ldB);

void LinearSolveHPD(int Size, int Nrhs, std::complex<double>* A,
                    int ldA, std::complex<double>* B, int ldB);

void LinearSolve(int Size, int Nrhs, double* A, int ldA, double* B, int ldB);

void LeastSquares(int Size, int Size2, int Nrhs, double* A, int ldA, double* B, int ldB);

void SingularValueDecomposition(int Size1, int Size2, double* A, double* U,
                                double* D, double* VT);

void SingularValueDecomposition(int Size1, int Size2,
                                std::complex<double>* A, std::complex<double>* U,
                                double* D, std::complex<double>* VH);

void SingularValueDecompositionFull(int Size1, int Size2, double* A, double* U, double* D, double* VT);

void SingularValueDecompositionFull(int Size1, int Size2,
                                    std::complex<double>* A, std::complex<double>* U,
                                    double* D, std::complex<double>* VH);

void EigenvaluesSymmetric(int Size, double* Data, int LeadingDim, double* Eigen);

void EigenvaluesHermitian(int Size, std::complex<double>* Data,
                          int LeadingDim, double* Eigen);

void EigenvaluesComplex(int Size, std::complex<double>* Data,
                        int LeadingDim, std::complex<double>* Eigen);

void DiagonalizeSymmetric(int Size, double* Data, int LeadingDim, double* Eigen);

void DiagonalizeHermitian(int Size, std::complex<double>* Data,
                          int LeadingDim, double* Eigen);

void Diagonalize(int Size, std::complex<double> const* Data,
                 int LeadingDim, std::complex<double>* Eigen,
                 std::complex<double>* LeftVectors, int LeftLeadingDim,
                 std::complex<double>* RightVectors, int RightLeadingDim);

void GeneralizedEigenSymmetric(int Size, double* A, int ldA, double* B, int ldB,
                               int First, int Last, double* Eigenval, double* Z,
                               int ldZ, double abstol);

void GeneralizedEigenHermitian(int Size, std::complex<double>* A, int ldA,
                               std::complex<double>* B, int ldB,
                               int First, int Last, double* Eigenval,
                               std::complex<double>* Z, int ldZ,
                               double abstol);

void TridiagonalizeHermitian(int Size, std::complex<double>* A, int ldA,
                             double* Diag,
                             double* SubDiag);

void CholeskyUpper(int Size, std::complex<double>* A, int ldA);
void CholeskyUpper(int Size, double* A, int ldA);

void CholeskyLower(int Size, std::complex<double>* A, int ldA);
void CholeskyLower(int Size, double* A, int ldA);

void InvertHPD(int Size, std::complex<double>* A, int ldA);

void InvertGeneral(int Size, std::complex<double>* A, int ldA);

void InvertUpperTriangular(int Size, std::complex<double>* A, int ldA);
void InvertLowerTriangular(int Size, std::complex<double>* A, int ldA);

void LQ_Factorize(int Size1, int Size2, double* A, int ldA, double* Tau);
void LQ_Factorize(int Size1, int Size2, std::complex<double>* A, int ldA, std::complex<double>* Tau);

void LQ_Construct(int Size1, int Size2, int k, double* A, int ldA, double* Tau);
void LQ_Construct(int Size1, int Size2, int k, std::complex<double>* A, int ldA, std::complex<double>* Tau);

} // namespace Private

// LinearSolveSPD

template <typename M1, typename M2, typename M1i, typename M2i>
struct ImplementLinearSolveSPD<M1, M2,
                               Concepts::MatrixExpression<double, M1i>,
                               Concepts::MatrixExpression<double,M2i>>
{
   typedef Matrix<double, ColMajor> result_type;
   static result_type apply(M1 const& m, M2 const& rhs)
   {
      PRECONDITION_EQUAL(size1(m), size2(m));
      PRECONDITION_EQUAL(size1(rhs), size1(m));

      Matrix<double> TempM(m);
      Matrix<double, ColMajor> TempRhs(rhs);

      // We don't care here whether A is in row-major or column major format,
      // since it is symmetric.
      CHECK(is_blas_matrix(TempM));
      // It does matter whether TempRhs is row- or column major, because we want the final vectors
      // to be column vectors.
      CHECK(is_blas_matrix(TempRhs));
      CHECK_EQUAL(stride1(TempRhs), 1);

      //   DEBUG_PRECONDITION(is_symmetric(m));
      //   DEBUG_PRECONDITION(min(EigenvaluesSymmetric(m)) > 0.0);

      Private::LinearSolveSPD(size1(TempM), size2(TempRhs), data(TempM), leading_dimension(TempM),
                              data(TempRhs), stride2(TempRhs) );
      return TempRhs;
   }
};

// LinearSolveHPD

template <typename M1, typename M2>
inline
void ImplementLinearSolveHPD(M1& m, M2& rhs)
{
   PRECONDITION_EQUAL(size1(m), size2(m));
   PRECONDITION_EQUAL(size1(rhs), size1(m));

   // m is hermitian, so it makes a difference whether it is row- or column-major.
   // we require column major here.
   PRECONDITION(is_blas_matrix(m));
   PRECONDITION_EQUAL(stride1(m), 1);
   // It does matter whether B is row- or column major, because we want the final vectors
   // to be column vectors.
   PRECONDITION(is_blas_matrix(rhs));
   PRECONDITION_EQUAL(stride1(rhs), 1);

   //   DEBUG_PRECONDITION(is_hermitian(m));
   //   DEBUG_PRECONDITION(min(EigenvaluesHermitian(m)) > 0.0);

   Private::LinearSolveHPD(size1(m), size2(rhs), data(m), leading_dimension(m),
                           data(rhs), stride2(rhs) );
}

template <typename M1, typename M2, typename T1, typename T2>
inline
Matrix<double, ColMajor>
LinearSolveHPD(M1 const& m, M2 const& rhs,
                    Concepts::MatrixExpression<double, T1>,
                    Concepts::MatrixExpression<double, T2>)
{
   Matrix<double, ColMajor> TempM(m);
   Matrix<double, ColMajor> TempRhs(rhs);
   return LinearSolveSPD(TempM, TempRhs);
}

template <typename M1, typename M2, typename T1, typename T2>
inline
Matrix<std::complex<double>, ColMajor>
LinearSolveHPD(M1 const& m, M2 const& rhs,
                    Concepts::MatrixExpression<std::complex<double>, T1>,
                    Concepts::MatrixExpression<std::complex<double>, T2>)
{
   Matrix<std::complex<double>, ColMajor> TempM(m);
   Matrix<std::complex<double>, ColMajor> TempRhs(rhs);
   ImplementLinearSolveHPD(TempM, TempRhs);
   return TempRhs;
}

template <typename M1, typename M2, typename T1, typename T2>
inline
Matrix<std::complex<double>, ColMajor>
LinearSolveHPD(M1 const& m, M2 const& rhs,
                    Concepts::MatrixExpression<double, T1>,
                    Concepts::MatrixExpression<std::complex<double>, T2>)
{
   Matrix<std::complex<double>, ColMajor> TempM(m);
   Matrix<std::complex<double>, ColMajor> TempRhs(rhs);
   ImplementLinearSolveHPD(TempM, TempRhs);
   return TempRhs;
}

template <typename M1, typename M2, typename T1, typename T2>
inline
Matrix<std::complex<double>, ColMajor>
LinearSolveHPD(M1 const& m, M2 const& rhs,
                    Concepts::MatrixExpression<std::complex<double>, T1>,
                    Concepts::MatrixExpression<double, T2>)
{
   Matrix<std::complex<double>, ColMajor> TempM(m);
   Matrix<std::complex<double>, ColMajor> TempRhs(rhs);
   ImplementLinearSolveHPD(TempM, TempRhs);
   return TempRhs;
}


template <typename M1, typename M2>
inline
Matrix<std::complex<double>, ColMajor>
LinearSolveHPD(M1 const& m, M2 const& rhs,
               typename boost::enable_if<is_matrix<M1> >::type*,
               typename boost::enable_if<is_matrix<M2> >::type*)
{
   return LinearSolveHPD(m, rhs,
                         typename interface<M1>::type(),
                         typename interface<M2>::type());
}

// LinearSolve

template <typename M1, typename M2>
inline
void ImplementLinearSolve(M1& m, M2& rhs)
{
   //PRECONDITION_EQUAL(size1(m), size2(m));
   PRECONDITION_EQUAL(size1(rhs), size1(m));

   // we require column major here.
   PRECONDITION(is_blas_matrix(m));
   PRECONDITION_EQUAL(stride1(m), 1);
   PRECONDITION(is_blas_matrix(rhs));
   PRECONDITION_EQUAL(stride1(rhs), 1);

   //   DEBUG_PRECONDITION(is_hermitian(m));
   //   DEBUG_PRECONDITION(min(EigenvaluesHermitian(m)) > 0.0);

   Private::LinearSolve(size1(m), size2(rhs), data(m), leading_dimension(m),
                        data(rhs), stride2(rhs) );
}

template <typename M1, typename M2, typename T1, typename T2>
inline
Matrix<double, ColMajor>
LinearSolve(M1 const& m, M2 const& rhs,
            Concepts::MatrixExpression<double, T1>,
            Concepts::MatrixExpression<double, T2>)
{
   Matrix<double, ColMajor> TempM(m);
   Matrix<double, ColMajor> TempRhs(rhs);
   ImplementLinearSolve(TempM, TempRhs);
   return TempRhs;
}

template <typename M1, typename M2>
inline
Matrix<double, ColMajor>
LinearSolve(M1 const& m, M2 const& rhs,
               typename boost::enable_if<is_matrix<M1> >::type*,
               typename boost::enable_if<is_matrix<M2> >::type*)
{
   return LinearSolve(m, rhs,
                      typename interface<M1>::type(),
                      typename interface<M2>::type());
}


template <typename M, typename V>
inline
void ImplementLeastSquares(M& m, V& rhs)
{
   //PRECONDITION_EQUAL(size1(m), size2(m));
   PRECONDITION_EQUAL(size(rhs), size1(m));

   // we require column major here.
   PRECONDITION(is_blas_matrix(m));
   PRECONDITION_EQUAL(stride1(m), 1);

   //   DEBUG_PRECONDITION(is_hermitian(m));
   //   DEBUG_PRECONDITION(min(EigenvaluesHermitian(m)) > 0.0);

   Private::LeastSquares(size1(m), size2(m), 1, data(m), leading_dimension(m),
                         data(rhs), size(rhs) );
}

template <typename M, typename V, typename T1, typename T2>
inline
Vector<double>
LeastSquares(M const& m, V const& rhs,
            Concepts::MatrixExpression<double, T1>,
            VECTOR_EXPRESSION(double, T2))
{
   Matrix<double, ColMajor> TempM(m);
   Vector<double> TempRhs(rhs);
   ImplementLeastSquares(TempM, TempRhs);
   Vector<double> Result(size2(m));
   Result = TempRhs[LinearAlgebra::range(0, size2(m))];
   return Result;
}

template <typename M, typename V>
inline
Vector<double>
LeastSquares(M const& m, V const& rhs,
               typename boost::enable_if<is_matrix<M> >::type*,
               typename boost::enable_if<is_vector<V> >::type*)
{
   return LeastSquares(m, rhs,
                       typename interface<M>::type(),
                       typename interface<V>::type());
}

// DiagonalizeSymmetric

template <typename M, typename Mi>
struct ImplementDiagonalizeSymmetric<M, Concepts::ContiguousMatrix<double, RowMajor, Mi>>
{
   typedef Vector<double> result_type;
   static result_type apply(M& m)
   {
      PRECONDITION_EQUAL(size1(m), size2(m));
      //      DEBUG_PRECONDITION(is_symmetric(M));        // this is an expensive test
      result_type Result(size1(m));
      DEBUG_CHECK_EQUAL(stride2(m),1);
      DEBUG_CHECK_EQUAL(difference_type(stride1(m)), difference_type(size2(m)));
      Private::DiagonalizeSymmetric(size1(m), data(m), stride1(m), data(Result));
      return Result;
   }
};

template <typename M, typename Mi>
struct ImplementDiagonalizeSymmetric<M, Concepts::StrideMatrix<double, RowMajor, Mi>>
{
   typedef Vector<double> result_type;
   static result_type apply(M& m)
   {
      if (int(stride2(m)) != 1 || int(stride1(m)) != int(size2(m)))
      {
         Matrix<double, RowMajor> Temp(m);
         result_type Result = DiagonalizeSymmetric(Temp);
         assign(m, Temp);
         return Result;
      }
      // else
      PRECONDITION_EQUAL(size1(m), size2(m));
      //      DEBUG_PRECONDITION(is_symmetric(M));        // this is an expensive test
      result_type Result(size1(m));
      DEBUG_CHECK_EQUAL(stride2(m),1);
      DEBUG_CHECK_EQUAL(int(stride1(m)),int(size2(m)));
      Private::DiagonalizeSymmetric(size1(m), data(m), stride1(m), data(Result));
      return Result;
   }
};

template <typename M, typename Mi>
struct ImplementDiagonalizeSymmetric<M, Concepts::MatrixExpression<double, Mi>>
{
   typedef Vector<double> result_type;
   static result_type apply(M& m)
   {
      Matrix<double, RowMajor> Temp(m);
      result_type Result = DiagonalizeSymmetric(Temp);
      assign(m, Temp);
      return Result;
   }
};

// DiagonalizeHermitian

template <typename M, typename Mi>
struct ImplementDiagonalizeHermitian<M, Concepts::MatrixExpression<double, Mi>>
   : ImplementDiagonalizeSymmetric<M> {};


template <typename M, typename Mi>
struct ImplementDiagonalizeHermitian<
   M
 ,Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>
>
{
   typedef Vector<double> result_type;
   static result_type apply(M& m)
   {
      PRECONDITION_EQUAL(size1(m), size2(m));
      //      DEBUG_PRECONDITION(is_hermitian(M));        // this is an expensive test
      result_type Result(size1(m));
      DEBUG_CHECK_EQUAL(stride2(m),1);
      DEBUG_CHECK_EQUAL(stride1(m), size2(m));
      Private::DiagonalizeHermitian(size1(m), data(m), stride1(m), data(Result));
      return Result;
   }
};

template <typename M, typename Mi>
struct ImplementDiagonalizeHermitian<
   M
 , Concepts::StrideMatrix<std::complex<double>, RowMajor, Mi>
>
{
   typedef Vector<double> result_type;
   static result_type apply(M& m)
   {
      if (difference_type(stride2(m)) != 1 || difference_type(stride1(m)) != difference_type(size2(m)))
      {
         Matrix<std::complex<double>, RowMajor> Temp(m);
         result_type Result = DiagonalizeHermitian(Temp);
         assign(m, Temp);
         return Result;
      }
      // else
      PRECONDITION_EQUAL(size1(m), size2(m));
      //      DEBUG_PRECONDITION(is_hermitian(M));        // this is an expensive test
      result_type Result(size1(m));
      DEBUG_CHECK_EQUAL(stride2(m),1);
      DEBUG_CHECK_EQUAL(std::size_t(stride1(m)),size2(m));
      Private::DiagonalizeHermitian(size1(m), data(m), stride1(m), data(Result));
      return Result;
   }
};

template <typename M, typename Mi>
struct ImplementDiagonalizeHermitian<M, Concepts::MatrixExpression<std::complex<double>, Mi>>
{
   typedef Vector<double> result_type;
   static result_type apply(M& m)
   {
      Matrix<std::complex<double>, RowMajor> Temp(m);
      result_type Result = DiagonalizeHermitian(Temp);
      assign(m, Temp);
      return Result;
   }
};

// Diagonalize

template <typename M, typename L, typename R, typename Mi, typename Li, typename Ri>
struct ImplementDiagonalize<M, L&, R&,
                            Concepts::MatrixExpression<std::complex<double>, Mi>,
                            Concepts::MatrixExpression<std::complex<double>, Li>,
                            Concepts::MatrixExpression<std::complex<double>, Ri>>
{
   typedef Vector<std::complex<double> > result_type;
   static result_type apply(M const& m, L& l, R& r)
   {
      PRECONDITION_EQUAL(size1(m), size2(m));

      std::size_t Size = size1(m);
      Matrix<std::complex<double>, RowMajor> Temp(m);
      Matrix<std::complex<double>, RowMajor> Left(Size, Size);
      Matrix<std::complex<double>, RowMajor> Right(Size, Size);
      Vector<std::complex<double> > Eigen(Size);

      Private::Diagonalize(Size, data(Temp), Size, data(Eigen),
                           data(Left), Size, data(Right), Size);

      r = conj(Left);    // interchange and conjugate, since we took the transpose matrix
      l = conj(Right);
      return Eigen;
   }
};

// SingularValueDecomposition

// TODO: we could avoid some temporaries here with some overloads for CONTIGUOUS_MATRIX
// But that is hardly necessary with the reference counted implementation

template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti>
struct ImplementSingularValueDecomposition<A, U, D, Vt,
                                           Concepts::MatrixExpression<double, Ai>,
                                           Concepts::MatrixExpression<double, Ui>,
                                           VECTOR_EXPRESSION(double, Di),
                                           Concepts::MatrixExpression<double, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d, Vt& vt) const
   {
      int m = a.size2();
      int n = a.size1();
      int min_mn = std::min(m,n);

      try_resize(u, n, min_mn);
      try_resize(d, min_mn);
      try_resize(vt, min_mn, m);

      if (min_mn == 0) return;

      Matrix<double> Acopy(a);
      Matrix<double> Ures(n, min_mn);
      Vector<double> Dres(min_mn);
      Matrix<double> Vtres(min_mn, m);

      Private::SingularValueDecomposition(size1(Acopy), size2(Acopy), data(Acopy),
                                          data(Ures), data(Dres), data(Vtres));
#if defined(RANDOMIZE_VECTORS)
      for (unsigned i = 0; i < min_mn; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Ures(LinearAlgebra::all,i) *= Phase;
         Vtres(i,LinearAlgebra::all) *= LinearAlgebra::conj(Phase);
      }
#endif
      assign(u, Ures);
      assign(d, Dres);
      assign(vt, Vtres);
   }
};

template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti>
struct ImplementSingularValueDecomposition<A, U, D, Vt,
                                           Concepts::MatrixExpression<std::complex<double>, Ai>,
                                           Concepts::MatrixExpression<std::complex<double>, Ui>,
                                           VECTOR_EXPRESSION(double, Di),
                                           Concepts::MatrixExpression<std::complex<double>, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d, Vt& vt) const
   {
      int m = a.size2();
      int n = a.size1();
      int min_mn = std::min(m,n);

      try_resize(u, n, min_mn);
      try_resize(d, min_mn);
      try_resize(vt, min_mn, m);

      if (min_mn == 0) return;

      Matrix<std::complex<double> > Acopy(a);
      Matrix<std::complex<double> > Ures(n, min_mn);
      Vector<double> Dres(min_mn);
      Matrix<std::complex<double> > Vtres(min_mn, m);

      Private::SingularValueDecomposition(size1(Acopy), size2(Acopy), data(Acopy),
                                          data(Ures), data(Dres), data(Vtres));
#if defined(RANDOMIZE_VECTORS)
      for (int i = 0; i < min_mn; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Ures(LinearAlgebra::all,i) *= Phase;
         Vtres(i,LinearAlgebra::all) *= LinearAlgebra::conj(Phase);
      }
#endif
      assign(u, Ures);
      assign(d, Dres);
      assign(vt, Vtres);
   }
};

template <typename A, typename U, typename D,
          typename Ai, typename Ui, typename Di>
struct ImplementSingularValueDecompositionLeft<A, U, D,
                                           Concepts::MatrixExpression<std::complex<double>, Ai>,
                                           Concepts::MatrixExpression<std::complex<double>, Ui>,
                                           VECTOR_EXPRESSION(double, Di)>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d) const
   {
      int m = a.size2();
      int n = a.size1();
      int min_mn = std::min(m,n);

      try_resize(u, n, min_mn);
      try_resize(d, min_mn);

      if (min_mn == 0) return;

      Matrix<std::complex<double> > Acopy(a);
      Matrix<std::complex<double> > Ures(n, min_mn);
      Vector<double> Dres(min_mn);

      Private::SingularValueDecomposition(size1(Acopy), size2(Acopy), data(Acopy), data(Ures), data(Dres), nullptr);
#if defined(RANDOMIZE_VECTORS)
      for (int i = 0; i < min_mn; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Ures(LinearAlgebra::all,i) *= Phase;
      }
#endif
      assign(u, Ures);
      assign(d, Dres);
   }
};

template <typename A, typename D, typename Vt,
          typename Ai, typename Di, typename Vti>
struct ImplementSingularValueDecompositionRight<A, D, Vt,
                                           Concepts::MatrixExpression<std::complex<double>, Ai>,
                                           VECTOR_EXPRESSION(double, Di),
                                           Concepts::MatrixExpression<std::complex<double>, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, D& d, Vt& vt) const
   {
      int m = a.size2();
      int n = a.size1();
      int min_mn = std::min(m,n);

      try_resize(d, min_mn);
      try_resize(vt, min_mn, m);

      if (min_mn == 0) return;

      Matrix<std::complex<double> > Acopy(a);
      Vector<double> Dres(min_mn);
      Matrix<std::complex<double> > Vtres(min_mn, m);

      Private::SingularValueDecomposition(size1(Acopy), size2(Acopy), data(Acopy), nullptr, data(Dres), data(Vtres));
#if defined(RANDOMIZE_VECTORS)
      for (int i = 0; i < min_mn; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Vtres(i,LinearAlgebra::all) *= LinearAlgebra::conj(Phase);
      }
#endif
      assign(d, Dres);
      assign(vt, Vtres);
   }
};


// version taking a diagonal matrix for D

// Real is either double or complex<double>
template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti, typename Real>
struct ImplementSingularValueDecomposition<A, U, D, Vt,
                                           Concepts::MatrixExpression<Real, Ai>,
                                           Concepts::MatrixExpression<Real, Ui>,
                                           Concepts::DiagonalMatrix<double, Di>,
                                           Concepts::MatrixExpression<Real, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d, Vt& vt) const
   {
      int m = a.size2();
      int n = a.size1();
      int min_mn = std::min(m,n);

      try_resize(u, n, min_mn);
      try_resize(d, min_mn, min_mn);
      try_resize(vt, min_mn, m);

      if (min_mn == 0) return;

      Matrix<Real> Acopy(a);
      Matrix<Real> Ures(n, min_mn);
      Vector<double> Dres(min_mn);
      Matrix<Real> Vtres(min_mn, m);

      Private::SingularValueDecomposition(size1(Acopy), size2(Acopy), data(Acopy),
                                          data(Ures), data(Dres), data(Vtres));
#if defined(RANDOMIZE_VECTORS)
      for (int i = 0; i < min_mn; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Ures(LinearAlgebra::all,i) *= Phase;
         Vtres(i,LinearAlgebra::all) *= LinearAlgebra::conj(Phase);
      }
#endif
      assign(u, Ures);
      assign(d.diagonal(), Dres);
      assign(vt, Vtres);
   }
};

// Real is either double or complex<double>
template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti, typename Real>
struct ImplementSingularValueDecompositionFull<A, U, D, Vt,
                                               Concepts::MatrixExpression<Real, Ai>,
                                               Concepts::MatrixExpression<Real, Ui>,
                                               Concepts::DiagonalMatrix<double, Di>,
                                               Concepts::MatrixExpression<Real, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d, Vt& vt) const
   {
      int m = a.size1();
      int n = a.size2();
      int min_mn = std::min(m,n);
      int max_mn = std::max(m,n);

      try_resize(u, m, max_mn);
      try_resize(d, max_mn, max_mn);
      try_resize(vt, max_mn, n);

      // If there are no non-zero singular values then we can't call LAPACK, and the usual SVD
      // isn't well defined.  But we can interpret the full SVD for an Mx0 matrix as an MxM random orthogonal
      // U, an MxM zero matrix D, and an Mx0 matrix Vt, and similarly for a 0xN matrix.
      if (min_mn == 0)
      {
         if (n == 0)
         {
            // set u to be a random m x m orthogonal matrix
            u = LinearAlgebra::random_unitary<Real>(m, m);
         }
         else // m == 0
         {
            // set vt to be a random n x n orthogonal matrix
            vt = LinearAlgebra::random_unitary<Real>(n, n);
         }
         zero_all(d);
         return;
      }

      Matrix<Real> Acopy(a);
      Matrix<Real> Ures(m, m);
      Vector<double> Dres(min_mn);
      Matrix<Real> Vtres(n, n);

      Private::SingularValueDecompositionFull(size1(Acopy), size2(Acopy), data(Acopy),
                                              data(Ures), data(Dres), data(Vtres));
      // randomize signs
#if defined(RANDOMIZE_VECTORS)
      for (int i = 0; i < min_mn; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Ures(LinearAlgebra::all,i) *= Phase;
         Vtres(i,LinearAlgebra::all) *= LinearAlgebra::conj(Phase);
      }
#endif
      zero_all(u);
      zero_all(d);
      zero_all(vt);
      assign(u(LinearAlgebra::all, LinearAlgebra::range(0, m)), Ures);
      assign(d.diagonal()[LinearAlgebra::range(0, min_mn)], Dres);
      assign(vt(LinearAlgebra::range(0, n), LinearAlgebra::all), Vtres);
   }
};

template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti, typename Real>
struct ImplementSingularValueDecompositionFullLeft<A, U, D, Vt,
                                                   Concepts::MatrixExpression<Real, Ai>,
                                                   Concepts::MatrixExpression<Real, Ui>,
                                                   Concepts::DiagonalMatrix<double, Di>,
                                                   Concepts::MatrixExpression<Real, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d, Vt& vt) const
   {
      if (a.size1() <= a.size2())
         return SingularValueDecomposition(a, u, d, vt);
      else
         return SingularValueDecompositionFull(a, u, d, vt);
   }
};

template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti, typename Real>
struct ImplementSingularValueDecompositionFullRight<A, U, D, Vt,
                                                   Concepts::MatrixExpression<Real, Ai>,
                                                   Concepts::MatrixExpression<Real, Ui>,
                                                   Concepts::DiagonalMatrix<double, Di>,
                                                   Concepts::MatrixExpression<Real, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d, Vt& vt) const
   {
      if (a.size2() <= a.size1())
         return SingularValueDecomposition(a, u, d, vt);
      else
         return SingularValueDecompositionFull(a, u, d, vt);
   }
};

template <typename A, typename U, typename D,
          typename Ai, typename Ui, typename Di>
struct ImplementSingularValueDecompositionLeftFull<A, U, D,
                                           Concepts::MatrixExpression<std::complex<double>, Ai>,
                                           Concepts::MatrixExpression<std::complex<double>, Ui>,
                                           VECTOR_EXPRESSION(double, Di)>
{
   typedef void result_type;
   void operator()(A const& a, U& u, D& d) const
   {
      int m = a.size1();
      int n = a.size2();

      try_resize(u, m, m);
      try_resize(d, m);

      // If there are no non-zero singular values then we can't call LAPACK, and the usual SVD
      // isn't well defined.  But we can interpret the full SVD for an Mx0 matrix as an MxM random orthogonal
      // U, an MxM zero matrix D, and an Mx0 matrix Vt, and similarly for a 0xN matrix.
      if (n == 0)
      {
         if (m > 0)
         {
            // set u to be a random m x m orthogonal matrix
            assign(u, LinearAlgebra::random_unitary<std::complex<double>>(m, m));
            zero_all(d);
         }
         //else we have a 0 x 0 matrix
         return;
      }

      Matrix<std::complex<double>> Acopy(a);
      Matrix<std::complex<double>> Ures(m, m);
      Vector<double> Dres(m);
      //TRACE(m)(n);

      if (m < n)
         Private::SingularValueDecomposition(size1(Acopy), size2(Acopy), data(Acopy), data(Ures), data(Dres), nullptr);
      else
         Private::SingularValueDecompositionFull(size1(Acopy), size2(Acopy), data(Acopy), data(Ures), data(Dres), nullptr);

#if defined(RANDOMIZE_VECTORS)
      for (int i = 0; i < m; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Ures(LinearAlgebra::all,i) *= Phase;
      }
#endif
      assign(u, Ures);
      assign(d, Dres);
   }
};

template <typename A, typename D, typename Vt,
          typename Ai, typename Di, typename Vti>
struct ImplementSingularValueDecompositionRightFull<A, D, Vt,
                                           Concepts::MatrixExpression<std::complex<double>, Ai>,
                                           VECTOR_EXPRESSION(double, Di),
                                           Concepts::MatrixExpression<std::complex<double>, Vti>>
{
   typedef void result_type;
   void operator()(A const& a, D& d, Vt& vt) const
   {
      int m = a.size1();
      int n = a.size2();

      try_resize(d, n);
      try_resize(vt, n, n);

      // If there are no non-zero singular values then we can't call LAPACK, and the usual SVD
      // isn't well defined.  But we can interpret the full SVD for an Mx0 matrix as an MxM random orthogonal
      // U, an MxM zero matrix D, and an Mx0 matrix Vt, and similarly for a 0xN matrix.
      if (m == 0)
      {
         if (n > 0)
         {
            // set vt to be a random n x n orthogonal matrix
            vt = LinearAlgebra::random_unitary<std::complex<double>>(n, n);
            zero_all(d);
         }
         return;
      }

      Matrix<std::complex<double> > Acopy(a);
      Vector<double> Dres(n);
      Matrix<std::complex<double> > Vtres(n, n);

      if (n < m)
         Private::SingularValueDecomposition(size1(Acopy), size2(Acopy), data(Acopy), nullptr, data(Dres), data(Vtres));
      else
         Private::SingularValueDecompositionFull(size1(Acopy), size2(Acopy), data(Acopy), nullptr, data(Dres), data(Vtres));

#if defined(RANDOMIZE_VECTORS)
      for (int i = 0; i < n; ++i)
      {
         double Phase = (rand() % 2) * 2.0 - 1.0;
         Vtres(i,LinearAlgebra::all) *= Phase;
      }
#endif
      assign(d, Dres);
      assign(vt, Vtres);
   }
};

// GeneralizedEigenSymmetric

template <typename A, typename B, typename Eigenval, typename Eigenvec,
          typename Ai, typename Bi, typename Eigenvali, typename Eigenveci>
struct ImplementGeneralizedEigenSymmetric<A, B, Eigenval, Eigenvec,
                                                             Concepts::MatrixExpression<double, Ai>,
                                                             Concepts::MatrixExpression<double, Bi>,
                                                             VECTOR_EXPRESSION(double, Eigenvali),
                                                             Concepts::MatrixExpression<double, Eigenveci>>
                                                             {
   typedef void result_type;
   void operator()(A const& a, B const& b, Eigenval& eigenval, Eigenvec& eigenvec,
                   Range const& Which, double abstol)
   {
      size_type N = size1(a);
      PRECONDITION_EQUAL(size2(a), N);
      PRECONDITION_EQUAL(size1(b), N);
      PRECONDITION_EQUAL(size2(b), N);
      PRECONDITION(Which.first() >= 0 && Which.last() <= N)(Which);
      DEBUG_PRECONDITION(is_symmetric(a))(a);
      DEBUG_PRECONDITION(is_symmetric(b))(b);
      DEBUG_PRECONDITION(min(EigenvaluesSymmetric(b)) > 0.0);  // B must be positive-definite

      try_resize(eigenval, size(Which));
      try_resize(eigenvec, size(Which), N);

      if (size(Which) == 0) return;   // corner case: no eigenvalues requested

      Matrix<double> Acopy(a);
      Matrix<double> Bcopy(b);

      Vector<double> EVal(size(Which));
      Matrix<double> EVec(size(Which), N);

      Private::GeneralizedEigenSymmetric(N, data(Acopy), N, data(Bcopy), N,
                                         Which.first()+1, Which.last(),
                                         data(EVal),
                                         data(EVec), N,
                                         abstol);

      assign(eigenval, EVal);
      assign(eigenvec, EVec);
   }
};

template <typename A, typename B, typename Eigenval, typename Eigenvec,
          typename Ai, typename Bi, typename Eigenvali, typename Eigenveci>
struct ImplementGeneralizedEigenSymmetric<A, B, Eigenval, Eigenvec,
                                          Concepts::MatrixExpression<std::complex<double>, Ai>,
                                          Concepts::MatrixExpression<std::complex<double>, Bi>,
                                          VECTOR_EXPRESSION(double, Eigenvali),
                                          Concepts::MatrixExpression<std::complex<double>, Eigenveci>>
{
   typedef void result_type;
   void operator()(A const& a, B const& b, Eigenval& eigenval, Eigenvec& eigenvec,
                   Range const& Which, double abstol)
   {
      size_type N = size1(a);
      PRECONDITION_EQUAL(size2(a), N);
      PRECONDITION_EQUAL(size1(b), N);
      PRECONDITION_EQUAL(size2(b), N);
      PRECONDITION(Which.first() >= 0 && Which.last() <= N)(Which);
      DEBUG_PRECONDITION(is_hermitian(a))(a);
      DEBUG_PRECONDITION(is_hermitian(b))(b);
      DEBUG_PRECONDITION(min(EigenvaluesHermitian(b)) > 0.0);  // B must be positive-definite

      try_resize(eigenval, size(Which));
      try_resize(eigenvec, size(Which), N);

      if (size(Which) == 0) return;   // corner case: no eigenvalues requested

      Matrix<std::complex<double> > Acopy(a);
      Matrix<std::complex<double> > Bcopy(b);

      Vector<double> EVal(size(Which));
      Matrix<std::complex<double> > EVec(size(Which), N);

      Private::GeneralizedEigenHermitian(N, data(Acopy), N, data(Bcopy), N,
                                         Which.first()+1, Which.last(),
                                         data(EVal),
                                         data(EVec), N,
                                         abstol);

      assign(eigenval, EVal);
      assign(eigenvec, EVec);
   }
};

//
// TridiagonalizeHermitian
//

template <typename M, typename Mi>
struct ImplementTridiagonalizeHermitian<M&, Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>>
{
   typedef Matrix<double> result_type;
   result_type operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Matrix<double> Result(size1(m), size2(m), 0.0);
      Vector<double> Diag(size1(m));
      Vector<double> Sub(size1(m));
      Private::TridiagonalizeHermitian(size1(m), data(m), size1(m), data(Diag), data(Sub));
      for (int i = 0; i < int(size1(m))-1; ++i)
      {
         Result(i,i) = Diag[i];
         Result(i,i+1) = Sub[i];
         Result(i+1,i) = conj(Sub[i]);
      }
      Result(size1(m)-1, size1(m)-1) = Diag[size1(m)-1];
      return Result;
   }
};

//
// CholeskyFactorize
//

template <typename M, typename Mi>
struct ImplementCholeskyFactorizeUpper<M&, Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::CholeskyLower(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementCholeskyFactorizeUpper<M&, Concepts::ContiguousMatrix<std::complex<double>, ColMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::CholeskyUpper(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementCholeskyFactorizeLower<M&, Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::CholeskyUpper(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementCholeskyFactorizeLower<M&, Concepts::ContiguousMatrix<std::complex<double>, ColMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::CholeskyLower(size1(m), data(m), size1(m));
   }
};

//
// Singular factorize
//

template <typename M, typename Mi>
struct ImplementSingularFactorize<M, Concepts::MatrixExpression<std::complex<double>, Mi>>
{
   typedef Matrix<std::complex<double> > result_type;
   result_type operator()(M const& m) const
   {
      result_type X(m);
      CHECK_EQUAL(size1(X), size2(X));
      Vector<double> Eigenvalues(size1(X));
      Private::DiagonalizeHermitian(size1(X), data(X), stride1(X), data(Eigenvalues));
      for (unsigned i = 0; i < size1(X); ++i)
      {
         X(i, all) *= std::sqrt(std::max(0.0, Eigenvalues[i]));
      }
      return X;
   }
};

template <typename M, typename Mi>
struct ImplementSingularFactorize<M, Concepts::MatrixExpression<double, Mi>>
{
   typedef Matrix<double> result_type;
   result_type operator()(M& m) const
   {
      result_type X(M);
      CHECK_EQUAL(size1(X), size2(X));
      Vector<double> Eigenvalues(size1(X));
      Private::DiagonalizeHermitian(size1(X), data(X), stride1(X), data(Eigenvalues));
      for (unsigned i = 0; i < size1(X); ++i)
      {
         X(i, all) *= std::sqrt(std::max(0.0, Eigenvalues[i]));
      }
      return X;
   }
};

//
// QR_Factorize
//

template <typename M, typename Mi>
struct ImplementQRFactorize<M&, Concepts::ContiguousMatrix<double, RowMajor, Mi>>
{
   typedef Matrix<double> result_type;
   result_type operator()(M& m) const
   {
      int s1 = size1(m);
      int s2 = size2(m);
      int sz = std::min(s1, s2);
      Vector<double> Tau(sz);
      Private::LQ_Factorize(size2(m), size1(m), data(m), stride1(m), data(Tau));

      // Convert the product of elementary reflectors into the Q matrix
      // NOTE: There is an alternative approach, which is to keep the matrix as it is
      // and use LAPACK dormqr() function to multiply another matrix directly, without explicitly constructing Q.
      // We could also use LAPACK dlacpy here
      Matrix<double> Q(s1, s1, 0.0);
      Q(LinearAlgebra::all, LinearAlgebra::range(0,sz)) = m(LinearAlgebra::range(0,s1), LinearAlgebra::range(0,sz));

      Private::LQ_Construct(s1, s1, sz, data(Q), stride1(Q), data(Tau));

      // Zero the unused parts of m, which now becomes upper-triangular
      for (int i = 0; i < s1; ++i)
      {
         int msz = std::min(i,s2);
         for (int j = 0; j < msz; ++j)
         {
            m(i,j) = 0.0;
         }
      }
      return Q;
   }
};

template <typename M, typename Mi>
struct ImplementQRFactorize<M&, Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>>
{
   typedef Matrix<std::complex<double>> result_type;
   result_type operator()(M& m) const
   {
      int s1 = size1(m);
      int s2 = size2(m);
      int sz = std::min(s1, s2);
      Vector<std::complex<double>> Tau(sz);
      Private::LQ_Factorize(size2(m), size1(m), data(m), stride1(m), data(Tau));

      // Convert the product of elementary reflectors into the Q matrix
      // NOTE: There is an alternative approach, which is to keep the matrix as it is
      // and use LAPACK dormqr() function to multiply another matrix directly, without explicitly constructing Q.
      // We could also use LAPACK dlacpy here
      Matrix<std::complex<double>> Q(s1, s1, 0.0);
      Q(LinearAlgebra::all, LinearAlgebra::range(0,sz)) = m(LinearAlgebra::range(0,s1), LinearAlgebra::range(0,sz));

      Private::LQ_Construct(s1, s1, sz, data(Q), stride1(Q), data(Tau));

      // Zero the unused parts of m, which now becomes upper-triangular
      for (int i = 0; i < s1; ++i)
      {
         int msz = std::min(i,s2);
         for (int j = 0; j < msz; ++j)
         {
            m(i,j) = 0.0;
         }
      }
      return Q;
   }
};

//
// InvertHPD
//

template <typename M, typename Orient, typename Mi>
struct ImplementInvertHPD<M&, Concepts::ContiguousMatrix<std::complex<double>, Orient, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      // We don't care here whether M is row- or column-major, it all sorts itself
      // out in the wash
      Private::InvertHPD(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertHPD<M&, Concepts::MatrixExpression<std::complex<double>, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      LinearAlgebra::Matrix<std::complex<double> > mc(m);
      InvertHPD(mc);
      m = mc;
   }
};

// TODO: we could implement a version for stride matrices...

//
// InvertGeneral
//

template <typename M, typename Orient, typename Mi>
struct ImplementInvertGeneral<M&, Concepts::ContiguousMatrix<std::complex<double>, Orient, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      // We don't care here whether M is row- or column-major, it all sorts itself
      // out in the wash
      Private::InvertGeneral(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertGeneral<M&, Concepts::MatrixExpression<std::complex<double>, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      LinearAlgebra::Matrix<std::complex<double> > mc(m);
      InvertGeneral(mc);
      m = mc;
   }
};

template <typename M>
Vector<double>
EigenvaluesSymmetric(M const& m)
{
   Matrix<double> Temp(m);
   return DiagonalizeSymmetric(Temp);
}

template <typename M>
Vector<double>
EigenvaluesHermitian(M const& m)
{
   Matrix<std::complex<double> > Temp(m);
   return DiagonalizeHermitian(Temp);
}

template <typename M>
Vector<std::complex<double>>
EigenvaluesComplex(M const& m)
{
   DEBUG_CHECK_EQUAL(size1(m), size2(m));
   Matrix<std::complex<double> > Temp(m);
   Vector<std::complex<double>> Result(size1(Temp));
   Private::EigenvaluesComplex(size1(Temp), data(Temp), size1(Temp), data(Result));
   return Result;
}

//
// Invert[Lower|Upper]Triangular
//

template <typename M, typename Mi>
struct ImplementInvertLowerTriangular<M&, Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertUpperTriangular(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertLowerTriangular<M&, Concepts::ContiguousMatrix<std::complex<double>, ColMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertLowerTriangular(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertUpperTriangular<M&, Concepts::ContiguousMatrix<std::complex<double>, RowMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertLowerTriangular(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertUpperTriangular<M&, Concepts::ContiguousMatrix<std::complex<double>, ColMajor, Mi>>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertUpperTriangular(size1(m), data(m), size1(m));
   }
};

} // namespace LinearAlgebra
