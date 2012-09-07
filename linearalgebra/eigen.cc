// -*- C++ -*- $Id$

#include "matrix.h"
#include "scalarmatrix.h"

namespace LinearAlgebra
{

namespace Private
{

void LinearSolveSPD(int Size, int Nrhs, double* A, int ldA, double* B, int ldB);

void LinearSolveHPD(int Size, int Nrhs, std::complex<double>* A, 
		    int ldA, std::complex<double>* B, int ldB);

void LinearSolve(int Size, int Nrhs, double* A, int ldA, double* B, int ldB);

void SingularValueDecomposition(int Size1, int Size2, double* A, double* U,
				double* D, double* VT);

void SingularValueDecomposition(int Size1, int Size2, 
				std::complex<double>* A, std::complex<double>* U,
				double* D, std::complex<double>* VH);

void EigenvaluesSymmetric(int Size, double* Data, int LeadingDim, double* Eigen);

void EigenvaluesHermitian(int Size, std::complex<double>* Data, 
			  int LeadingDim, double* Eigen);

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

void LQ_Factorize(int Size1, int Size2, std::complex<double>* A, int ldA, std::complex<double>* Tau);

} // namespace Private

// LinearSolveSPD

template <typename M1, typename M2, typename M1i, typename M2i>
struct ImplementLinearSolveSPD<M1, M2, 
			       MATRIX_EXPRESSION(double, M1i), 
			       MATRIX_EXPRESSION(double,M2i)>
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
		    MATRIX_EXPRESSION(double, T1), 
		    MATRIX_EXPRESSION(double, T2))
{
   Matrix<double, ColMajor> TempM(m);
   Matrix<double, ColMajor> TempRhs(rhs);
   return LinearSolveSPD(TempM, TempRhs);
}

template <typename M1, typename M2, typename T1, typename T2>
inline
Matrix<std::complex<double>, ColMajor> 
LinearSolveHPD(M1 const& m, M2 const& rhs, 
		    MATRIX_EXPRESSION(std::complex<double>, T1), 
		    MATRIX_EXPRESSION(std::complex<double>, T2))
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
		    MATRIX_EXPRESSION(double, T1), 
		    MATRIX_EXPRESSION(std::complex<double>, T2))
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
		    MATRIX_EXPRESSION(std::complex<double>, T1), 
		    MATRIX_EXPRESSION(double, T2))
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
   PRECONDITION_EQUAL(size1(m), size2(m));
   PRECONDITION_EQUAL(size1(rhs), size1(m));

   // m is hermitian, so it makes a difference whether it is row- or column-major.
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
            MATRIX_EXPRESSION(double, T1), 
            MATRIX_EXPRESSION(double, T2))
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



// DiagonalizeSymmetric

template <typename M, typename Mi>
struct ImplementDiagonalizeSymmetric<M, STRIDE_MATRIX(double, RowMajor, Mi)>
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
struct ImplementDiagonalizeSymmetric<M, CONTIGUOUS_MATRIX(double, RowMajor, Mi)>
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
struct ImplementDiagonalizeSymmetric<M, MATRIX_EXPRESSION(double, Mi)>
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
struct ImplementDiagonalizeHermitian<M, MATRIX_EXPRESSION(double, Mi)>
   : ImplementDiagonalizeSymmetric<M> {};


template <typename M, typename Mi>
struct ImplementDiagonalizeHermitian<
   M
 , STRIDE_MATRIX(std::complex<double>, RowMajor, Mi)
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
struct ImplementDiagonalizeHermitian<
   M
 ,CONTIGUOUS_MATRIX(std::complex<double>, RowMajor, Mi)
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
struct ImplementDiagonalizeHermitian<M, MATRIX_EXPRESSION(std::complex<double>, Mi)>
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
                            MATRIX_EXPRESSION(std::complex<double>, Mi),
                            MATRIX_EXPRESSION(std::complex<double>, Li),
                            MATRIX_EXPRESSION(std::complex<double>, Ri)>
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
                                           MATRIX_EXPRESSION(double, Ai),
                                           MATRIX_EXPRESSION(double, Ui),
                                           VECTOR_EXPRESSION(double, Di),
                                           MATRIX_EXPRESSION(double, Vti)>
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
      assign(u, Ures);
      assign(d, Dres);
      assign(vt, Vtres);
   }
};

template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti>
struct ImplementSingularValueDecomposition<A, U, D, Vt,
                                           MATRIX_EXPRESSION(std::complex<double>, Ai),
                                           MATRIX_EXPRESSION(std::complex<double>, Ui),
                                           VECTOR_EXPRESSION(double, Di),
                                           MATRIX_EXPRESSION(std::complex<double>, Vti)>
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
      assign(u, Ures);
      assign(d, Dres);
      assign(vt, Vtres);
   }
};

// version taking a diagonal matrix for D

template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti>
struct ImplementSingularValueDecomposition<A, U, D, Vt,
                                           MATRIX_EXPRESSION(double, Ai),
                                           MATRIX_EXPRESSION(double, Ui),
                                           DIAGONAL_MATRIX(double, Di),
                                           MATRIX_EXPRESSION(double, Vti)>
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
      assign(u, Ures);
      assign(d.diagonal(), Dres);
      assign(vt, Vtres);
   }
};

template <typename A, typename U, typename D, typename Vt,
          typename Ai, typename Ui, typename Di, typename Vti>
struct ImplementSingularValueDecomposition<A, U, D, Vt,
                                           MATRIX_EXPRESSION(std::complex<double>, Ai),
                                           MATRIX_EXPRESSION(std::complex<double>, Ui),
                                           DIAGONAL_MATRIX(double, Di),
                                           MATRIX_EXPRESSION(std::complex<double>, Vti)>
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
      assign(u, Ures);
      assign(d.diagonal(), Dres);
      assign(vt, Vtres);
   }
};

// GeneralizedEigenSymmetric

template <typename A, typename B, typename Eigenval, typename Eigenvec,
          typename Ai, typename Bi, typename Eigenvali, typename Eigenveci>
struct ImplementGeneralizedEigenSymmetric<A, B, Eigenval, Eigenvec,
                                          MATRIX_EXPRESSION(double, Ai),
                                          MATRIX_EXPRESSION(double, Bi),
                                          VECTOR_EXPRESSION(double, Eigenvali),
                                          MATRIX_EXPRESSION(double, Eigenveci)>
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
                                          MATRIX_EXPRESSION(std::complex<double>, Ai),
                                          MATRIX_EXPRESSION(std::complex<double>, Bi),
                                          VECTOR_EXPRESSION(double, Eigenvali),
                                          MATRIX_EXPRESSION(std::complex<double>, Eigenveci)>
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
struct ImplementTridiagonalizeHermitian<M&, CONTIGUOUS_MATRIX(std::complex<double>, RowMajor, Mi)>
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
struct ImplementCholeskyFactorizeUpper<M&, CONTIGUOUS_MATRIX(std::complex<double>, RowMajor, Mi)>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::CholeskyLower(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementCholeskyFactorizeUpper<M&, CONTIGUOUS_MATRIX(std::complex<double>, ColMajor, Mi)>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::CholeskyUpper(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementCholeskyFactorizeLower<M&, CONTIGUOUS_MATRIX(std::complex<double>, RowMajor, Mi)>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::CholeskyUpper(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementCholeskyFactorizeLower<M&, CONTIGUOUS_MATRIX(std::complex<double>, ColMajor, Mi)>
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
struct ImplementSingularFactorize<M, MATRIX_EXPRESSION(std::complex<double>, Mi)>
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
struct ImplementSingularFactorize<M, MATRIX_EXPRESSION(double, Mi)>
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
struct ImplementQRFactorize<M&, MATRIX_EXPRESSION(std::complex<double>, Mi)>
{
   typedef Matrix<std::complex<double> > result_type;
   result_type operator()(M const& m) const
   {
      Matrix<std::complex<double> > Mat(m);
      return QR_Factorize(m);
   }
};

template <typename M, typename Mi>
struct ImplementQRFactorize<M&, CONTIGUOUS_MATRIX(std::complex<double>, RowMajor, Mi)>
{
   typedef Matrix<std::complex<double> > result_type;
   result_type operator()(M& m) const
   {
      int s1 = size1(m);
      int s2 = size2(m);
      int sz = std::min(s1, s2);
      Vector<std::complex<double> > Tau(sz);
      Private::LQ_Factorize(size2(m), size1(m), data(m), stride1(m), data(Tau));

      Matrix<std::complex<double> > U = ScalarMatrix<double>(size1(m), size1(m), 1.0);
      Matrix<std::complex<double> > v(size1(m),1, 0.0);
      for (int i = sz-1; i >= 0; --i)
      {
         v(i,0) = 1.0;
         for (int j = i+1; j < s1; ++j)
         {
            v(j,0) = m(j,i); m(j,i) = 0.0;  // zero out the lower part while forming the elementary reflector
         } 
         U = Matrix<std::complex<double> >(ScalarMatrix<double>(size1(m), size1(m), 1.0) - conj(Tau[i])*v*herm(v)) * U;
      }
      return U;
   }
};

//
// InvertHPD
//

template <typename M, typename Orient, typename Mi>
struct ImplementInvertHPD<M&, CONTIGUOUS_MATRIX(std::complex<double>, Orient, Mi)>
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
struct ImplementInvertHPD<M&, MATRIX_EXPRESSION(std::complex<double>, Mi)>
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
struct ImplementInvertGeneral<M&, CONTIGUOUS_MATRIX(std::complex<double>, Orient, Mi)>
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
struct ImplementInvertGeneral<M&, MATRIX_EXPRESSION(std::complex<double>, Mi)>
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

//
// Invert[Lower|Upper]Triangular
//

template <typename M, typename Mi>
struct ImplementInvertLowerTriangular<M&, CONTIGUOUS_MATRIX(std::complex<double>, RowMajor, Mi)>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertUpperTriangular(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertLowerTriangular<M&, CONTIGUOUS_MATRIX(std::complex<double>, ColMajor, Mi)>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertLowerTriangular(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertUpperTriangular<M&, CONTIGUOUS_MATRIX(std::complex<double>, RowMajor, Mi)>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertLowerTriangular(size1(m), data(m), size1(m));
   }
};

template <typename M, typename Mi>
struct ImplementInvertUpperTriangular<M&, CONTIGUOUS_MATRIX(std::complex<double>, ColMajor, Mi)>
{
   typedef void result_type;
   void operator()(M& m) const
   {
      CHECK_EQUAL(size1(m), size2(m));
      Private::InvertUpperTriangular(size1(m), data(m), size1(m));
   }
};

} // namespace LinearAlgebra
