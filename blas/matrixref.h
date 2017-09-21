// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_BLAS_MATRIXREF_H)
#define MPTOOLKIT_BLAS_MATRIXREF_H

#include <list>
#include <mutex>
#include <iostream>

namespace blas
{

//
// number_traits
// basic information about scalar types that can be put in dense matrices
//

template <typename T>
struct number_traits
{
   using type = T;

   T zdero() constexpr { return 0; }
   T identity() constexpr { return 1; }
};

template <>
struct number_traits<float>
{
   using type = float;

   float zero() constexpr { return 0.0f; }
   float identity() constexpr { return 1.0f; }
};

template <>
struct number_traits<double>
{
   using type = double;

   double zero() constexpr { return 0.0; }
   double identity() constexpr { return 1.0; }
};

template <>
struct number_traits<std::complex<float>>
{
   using type = std::complex<float>;

   std::complex<float> zero() constexpr { return {0.0f,0.0f}; }
   std::complex<float> identity() constexpr { return {1.0f,0.0f}; }
};

template <>
struct number_traits<std::complex<double>>
{
   using type = std::complex<double>;

   std::complex<double> zero() constexpr { return {0.0,0.0}; }
   std::complex<double> identity() constexpr { return {1.0,0.0}; }
};

#if defined(HAVE_FLOAT128)
template <>
struct number_traits<float128>
{
   using type = float128;

   float128 zero() constexpr { return 0.0Q; }
   float128 identity() constexpr { return 1.0Q; }
};

template <>
struct number_traits<std::complex<float128>>
{
   using type = std::complex<float128>;

   std::complex<float128> zero() constexpr { return {0.0Q,0.0Q}; }
   std::complex<float128> identity() constexpr { return {1.0Q,0.0Q}; }
};
#endif

//
// MatrixRef : generic base class for a matrix using expression templates.
// MatrixRef is a reference to a matrix of type ValueType.  BaseType is 
// a concrete type, eg the two currently supported types are
// Matrix<ValueType> or gpu_matrix<ValueType>.
//
// If DerivedType == BaseType, then the matrix is actually a concrete matrix
// and can be constructed, used as a temporary, etc.
//

template <typename ValueType, BaseType, typename DerivedType = BaseType>
class MatrixRef<ValueType, BaseType, DerivedType>
{
   public:
      using value_type                        = ValueType;
      using derived_type                      = Derived;
      using base_type                         = Base;

      // default construction and move construction are defined, no copying.

      MatrixRef() = default;
      ~MatrixRef() = default;
      MatrixRef(MatrixRef&& Other) = default;
      MatrixRef(MatrixRef const&) = delete;
      MatrixRef& operator=(MatrixRef&&) = delete;
      MatrixRef& operator=(MatrixRef const&) = delete;

      derived_type& as_derived() { return *static_cast<derived_type*>(this); }
      derived_type const& as_derived() { return *static_cast<derived_type const*>(this); }

      int rows() const { return this->as_derived().rows(); }
      int cols() const { return this->as_derived().cols(); }
};

// Specialization of a MatrixRef for a matrix that can be directly used in BLAS-like calls,
// ie it physically exists, and has a valid BLAS trans parameter (ie, N, T, C, R).
// We expose the leading_dimension() function, but not the underlying pointer - that needs to be
// obtained via the derived class.

template <typename ValueType, typename BaseType, DerivedType = BaseType>
class BlasMatrix : public MatrixRef<ValueType, BaseType, DerivedType>
{
   public:
      using value_type   = ValueType;
      using derived_type = Derived;
      using base_type    = Base;

      MatrixRef() = default;
      ~MatrixRef() = default;
      MatrixRef(BlasMatrix&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }
      
      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation, )
      char trans() const { return this->as_derived().trans(); }
};

// proxy class for the transpose of a BlasMatrix

template <typename T, typename BaseType>
class BlasMatrixTrans : public BlasMatrix<T, BaseType, BlasMatrixTrans<T, BaseType>>
{
   public:
      using base_type = BaseType;
      using value_type    = T;

      BlasMatrixTrans(base_type const& Base_) : Base(Base_) {}

      int rows() const { return U.cols(); }
      int cols() const { return U.rows(); }
      int leading_dimension() const { return U.leading_dimension(); }
      char trans() const { return blas_traits<value_type>::trans_trans(U.trans()); }

      base_type const& base() const { return Base; }

   private:
      base_type const& Base;
};

// proxy class for the hermitian conjugate of a BlasMatrix

template <typename T, typename BaseType>
class BlasMatrixHerm : public BlasMatrix<T, BaseType, BlasMatrixHerm<T, BaseType>>
{
   public:
      using base_type = BaseType;
      using value_type    = T;

      BlasMatrixTrans(base_type const& Derived_) : Derived(D) {}

      int rows() const { return U.cols(); }
      int cols() const { return U.rows(); }
      int leading_dimension() const { return U.leading_dimension(); }
      char trans() const { return blas_traits<value_type>::trans_herm(U.trans()); }

      base_type const& base() const { return U; }

   private:
      base_type const& U;
};

// proxy class for the complex conjugate of a BlasMatrix

template <typename T, typename BaseType>
class BlasMatrixConj : public BlasMatrix<T, BaseType, BlasMatrixConj<T, BaseType>>
{
   public:
      using base_type = BaseType;
      using value_type    = T;

      BlasMatrixConj(base_type const& Derived_) : Derived(D) {}

      int rows() const { return U.rows(); }
      int cols() const { return U.cols(); }
      int leading_dimension() const { return U.leading_dimension(); }
      char trans() const { return blas_traits<value_type>::trans_conj(U.trans()); }

      base_type const& base() const { return U; }

   private:
      base_type const& U;
};

// trans

template <typename T, typename Base, typename Derived>
BlasMatrixTrans<T, Derived>
trans(BlasMatrix<T, Base, Derived> const& x)
{
   return BlasMatrixTrans<T, Derived>(x.as_derived());
}

template <typename T, typename Derived>
Dericed const&
trans(BlasMatrixTrans<T, Derived> const& x)
{
   return x.base();
}

template <typename T, typename U>
BlasMatrixConj<T, U>
trans(BlasMatrixHerm<T, U> const& x)
{
   return BlasMatrixConj<T, U>(x.base());
}

template <typename T, typename U>
BlasMatrixHerm<T, U>
trans(BlasMatrixConj<T, U> const& x)
{
   return BlasMatrixHerm<T, U>(x.base());
}

// herm

template <typename T, typename Base, typename Derived>
BlasMatrixHerm<T, Derived>
herm(BlasMatrix<T, Base, Derived> const& x)
{
   return BlasMatrixHerm<T, Derived>(x.as_derived());
}


template <typename T>
BlasMatrixHerm<T, U>
conj(BlasMatrixTrans<T, U> const& x)
{
   return BlasMatrixHerm<T, U>(x.base());
}

template <typename T, typename U>
U const&
herm(BlasMatrixHerm<T, U> const& x)
{
   return x.base();
}

template <typename T, typename U>
BlasMatrixTrans<T, U>
herm(BlasMatrixConj<T, U> const& x)
{
   return BlasMatrixTrans<T, U>(x.base());
}

template <typename T, typename U>
BlasMatrixConj<T, U>
herm(BlasMatrixTrans<T, U> const& x)
{
   return BlasMatrixConj<T, U>(x.base());
}

// specialization for real types, herm() reduces to trans()
template <typename Base, Derived>
BlasMatrixTrans<float, Derived>
herm(BlasMatrix<float, Base, Derived> const& x)
{
   return BlasMatrixHerm<float, Derived>(x.as_derived());
}

template <typename Base, Derived>
BlasMatrixTrans<double, Derived>
herm(BlasMatrix<double, Derived> const& x)
{
   return BlasMatrixHerm<double, U>(x.as_derived());
}

#if defined(HAVE_FLOAT128)
template <typename Base, Derived>
BlasMatrixTrans<float128, Derived>
herm(BlasMatrix<float128, Derived> const& x)
{
   return BlasMatrixHerm<float128, U>(x.as_derived());
}
#endif

// conj

template <typename T>
BlasMatrixConj<Matrix<T>, U>
conj(BlasMatrix<Matrix<T>, U> const& x)
{
   return BlasMatrixConj<Matrix<T>, U>(x.as_derived());
}

// specialization for real types, conj() reduces to a no-op
template <typename Base, Derived>
Derived&
conj(BlasMatrix<float, Base, Derived> const& x)
{
   return x.as_derived();
}

template <typename Base, Derived>
Derived&
conj(BlasMatrix<double, Base, Derived> const& x)
{
   return x.as_derived();
}

#if defined(HAVE_FLOAT128)
template <typename Base, Derived>
Derived&
conj(BlasMatrix<float128, Base, Derived> const& x)
{
   return x.as_derived();
}
#endif


// expression templates for gemm operation: alpha * op(A) * op(B)
// where op is normal, transpose, or hermitian, or conjugate

// firstly, for a product of a matrix and a scalar

template <typename T, typename U>
class ScalarMatrixProduct : public MatrixRef<T, U, ScalarMatrixProduct<T,U>>
{
   public:
      using base_type = U;

      ScalarMatrixProduct(T const& Factor_, U const& A) : Factor(Factor_), A(A_) {}

      int rows() const { return A.rows(); }
      int cols() const { return A.cols(); }

      base_type const& base() const { return A; }
      T factor() const { return Factor; }
   
   private:
      T Factor;
      U const& A;
};

template <typename T, typename U, typename Derived>
ScalarMatrixProduct<T, Derived>
operator*(T const& alpha, MatrixRef<T, U, Derived> const& M)
{
   return ScalarMatrixProduct<T, Derived>(alpha, M.as_derived());
}

void add(MatrixRef<T, BaseType>& A, ScalarMatrixProduct<T, U> const& B)
{
   add_scaled(A.as_derived(), B.factor(), B.base());
}

void subtract(MatrixRef<T, BaseType>& A, ScalarMatrixProduct<T, U> const& B)
{
   subtract_scaled(A.as_derived(), B.factor(), B.base());
}

template <typename T, typename BaseType, typename U, typename V>
struct MatrixProduct : public MatrixRef<T, BaseType, MatrixProduct<T, U, V>>
{
   MatrixProduct(MatrixRef<T, BaseType, V> const& A_, MatrixRef<T, W, X> cosnt& B_) : number_traits<T>::identity(), A(A_), B(B_) {}

   MatrixProduct(T const& Factor_, MatrixRef<T, U, V> const& A_, MatrixRef<T, W, X> cosnt& B_) 
      : Factor(Factor_), A(A_), B(B_) {}

   int rows() const { return A.rows(); }
   int cols() const { return B.cols(); }

   T Factor;
   AType const& A;
   BType const& B;
};

template <typename T, typename V, typename X>
MatrixProduct<T, U, V>
operator*(MatrixRef<T, U, V> const& A, MatrixRef<T, W, X> const& B)
{
   return MatrixProduct<T, U, V>(A.as_derived(), B.as_derived());
}



// BLAS gemm wrapper
template <typename T, typename U, typename V>
inline
void gemm(T alpha, BlasMatrix<Matrix<T>, U> const& A, BlasMatrix<Matrix<T>, V> const& B, T beta, Matrix<T>& C)
{
   gemm(A.trans(), B.trans(), alpha, A.rows(), A.cols(), B.cols(), alpha, A.as_derived().data(), A.leading_dim(), 
	B.as_derived().data(), B.leading_dim(), beta, C.data(), C.leading_dim());
}

// matrix product with a scalar, alpha * A * B

template <typename T, typename W, typename U, typename V>
void assign(MatrixRef<T>& C, MatrixRef<W, ScalarMatrixProduct<U, V> const& a)
{
   gemm(a.Operation.Factor, a.Operation.A, a.Operation.B, 0.0, C.as_derived());
}

template <typename T, typename W, typename U, typename V>
void add(MatrixRef<T>& C, MatrixRef<W, ScalarMatrixProduct<U, V> const& a)
{
   gemm(a.Operation.Factor, a.Operation.A, a.Operation.B, 1.0, C.as_derived());
}

template <typename T, typename U, typename V>
void subtract(Matrix<T>& C, MatrixRef<Matrix<T>, ScalarMatrixProduct<U, V>> const& a)
{
   gemm(a.Operation.Factor, a.Operation.A, a.Operation.B, -1.0, C);
}


template <typename T>
template <typename W, typename U>
Matrix<T>::Matrix(MatrixRef<W, U> const& E, arena Arene_)
   : Rows(E.rows()), Cols(E.cols()), LeadingDImension(leading_dimension(E.rows(), E.cols())),
     Arena(Arena_),
     Data(Arena.allocate(LeadingDimension * Cols, sizeof(T)))
{
   assign(*this, E.as_derived());
}


template <typename T>
void gemm(T alpha, Matrix<T> const& A, Matrix<T> const& B, T beta, Matrix<T>& C)
{
   DEBUC_CHECK_EQUAL(A.cols(), B.rows());
   gemm(operation::normal, operation::normal, A.rows(), A.cols(), B.cols(), alpha, A.data(), 
	A.leading_dim(), B.data(), B.leading_dim(), beta, C.data(), C.leading_dim());
}




} // namespace blas

#endif
