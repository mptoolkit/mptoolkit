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
#include "safe-conversions.h"
#include "number_traits.h"

namespace blas
{

//
// MatrixRef : generic base class for a matrix using expression templates.
// MatrixRef is a reference to a matrix of type ValueType.  BaseType is
// a concrete type, eg the two currently supported types are
// Matrix<ValueType> or gpu_matrix<ValueType>.
//
// If DerivedType == BaseType, then the matrix is actually a concrete matrix
// and can be constructed, used as a temporary, etc.
//

template <typename ValueType, typename BaseType, typename DerivedType = BaseType>
class MatrixRef
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using base_type      = BaseType;

      // default construction and move construction are defined, no copying.

      MatrixRef() = default;
      ~MatrixRef() = default;
      MatrixRef(MatrixRef&& Other) = default;
      MatrixRef(MatrixRef const&) = delete;
      MatrixRef& operator=(MatrixRef&&) = delete;
      MatrixRef& operator=(MatrixRef const&) = delete;

      derived_type& as_derived() { return *static_cast<derived_type*>(this); }
      derived_type const& as_derived() const { return *static_cast<derived_type const*>(this); }

      int rows() const { return this->as_derived().rows(); }
      int cols() const { return this->as_derived().cols(); }
};

// Specialization of a MatrixRef for a matrix that can be directly used in BLAS-like calls,
// ie it physically exists and can be addressed (although not necessarily in main memory,
// eg it might be on some other device such as a GPU), and has a valid BLAS trans parameter
// (ie, N, T, C, R).
// The BaseType is the concrete type, which allows us to specialize blas calls for different
// devices.
// We expose the leading_dimension() function, but not the memory storage, that needs
// to be done separately for each BaseType.

template <typename BaseTyepe>
struct blas_traits {};

template <typename ValueType, typename BaseType, typename DerivedType = BaseType>
class BlasMatrix : public MatrixRef<ValueType, BaseType, DerivedType>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using base_type          = BaseType;
      using storage_type       = typename blas_traits<base_type>::storage_type;
      using const_storage_type = typename blas_traits<base_type>::const_storage_type;

      BlasMatrix() = default;
      ~BlasMatrix() = default;
      BlasMatrix(BlasMatrix&& Other) = default;

      int leading_dimension() const { return this->as_derived().leading_dimension(); }

      // The BLAS trans parameter ('N', 'T', 'C', 'R' - note 'R' is a BLAS extension denoting conjugation)
      char trans() const { return this->as_derived().trans(); }

      storage_type storage() { return this->as_derived().storage(); }
      const_storage_type storage() const { return this->as_derived().storage(); }
};

// proxy class for the transpose of a BlasMatrix

template <typename T, typename BaseType>
class BlasMatrixTrans : public BlasMatrix<T, BaseType, BlasMatrixTrans<T, BaseType>>
{
   public:
      using base_type = BaseType;
      using value_type    = T;
      using storage_type       = typename blas_traits<base_type>::storage_type;
      using const_storage_type = typename blas_traits<base_type>::const_storage_type;

      BlasMatrixTrans(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.cols(); }
      int cols() const { return Base.rows(); }
      int leading_dimension() const { return Base.leading_dimension(); }
      char trans() const { return number_traits<value_type>::blas_trans(Base.trans()); }

      base_type const& base() const { return Base; }

      storage_type storage() { return Base.storage(); }
      const_storage_type storage() const { return Base.storage(); }

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
      using storage_type       = typename blas_traits<base_type>::storage_type;
      using const_storage_type = typename blas_traits<base_type>::const_storage_type;

      BlasMatrixHerm(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.cols(); }
      int cols() const { return Base.rows(); }
      int leading_dimension() const { return Base.leading_dimension(); }
      char trans() const { return number_traits<value_type>::blas_herm(Base.trans()); }

      base_type const& base() const { return Base; }

      storage_type storage() { return Base.storage(); }
      const_storage_type storage() const { return Base.storage(); }

   private:
      base_type const& Base;
};

// proxy class for the complex conjugate of a BlasMatrix

template <typename T, typename BaseType>
class BlasMatrixConj : public BlasMatrix<T, BaseType, BlasMatrixConj<T, BaseType>>
{
   public:
      using base_type = BaseType;
      using value_type    = T;
      using storage_type       = typename blas_traits<base_type>::storage_type;
      using const_storage_type = typename blas_traits<base_type>::const_storage_type;

      BlasMatrixConj(base_type const& Base_) : Base(Base_) {}

      int rows() const { return Base.rows(); }
      int cols() const { return Base.cols(); }
      int leading_dimension() const { return Base.leading_dimension(); }
      char trans() const { return number_traits<value_type>::blas_conj(Base.trans()); }
      T* data() { return Base.data(); }
      T const* data() const { return Base.data(); }

      base_type const& base() const { return Base; }

      storage_type storage() { return Base.storage(); }
      const_storage_type storage() const { return Base.storage(); }

   private:
      base_type const& Base;
};

// trans

template <typename T, typename Base, typename Derived>
BlasMatrixTrans<T, Derived>
trans(BlasMatrix<T, Base, Derived> const& x)
{
   return BlasMatrixTrans<T, Derived>(x.as_derived());
}

template <typename T, typename Derived>
Derived const&
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


template <typename T, typename U>
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
template <typename Base, typename Derived>
BlasMatrixTrans<float, Derived>
herm(BlasMatrix<float, Base, Derived> const& x)
{
   return BlasMatrixTrans<float, Derived>(x.as_derived());
}

template <typename Derived>
BlasMatrixTrans<double, Derived>
herm(BlasMatrix<double, Derived> const& x)
{
   return BlasMatrixTrans<double, Derived>(x.as_derived());
}

#if defined(HAVE_FLOAT128)
template <typename Base, typename Derived>
BlasMatrixTrans<float128, Derived>
herm(BlasMatrix<float128, Base, Derived> const& x)
{
   return BlasMatrixTrans<float128, Derived>(x.as_derived());
}
#endif

// conj

template <typename T, typename Base, typename Derived>
BlasMatrixConj<T, Derived>
conj(BlasMatrix<T, Base, Derived> const& x)
{
   return BlasMatrixConj<T, Derived>(x.as_derived());
}

// specialization for real types, conj() reduces to a no-op
template <typename Base, typename Derived>
Derived&
conj(BlasMatrix<float, Base, Derived> const& x)
{
   return x.as_derived();
}

template <typename Base, typename Derived>
Derived&
conj(BlasMatrix<double, Base, Derived> const& x)
{
   return x.as_derived();
}

#if defined(HAVE_FLOAT128)
template <typename Base, typename Derived>
Derived&
conj(BlasMatrix<float128, Base, Derived> const& x)
{
   return x.as_derived();
}
#endif

// assignment

template <typename T, typename BaseType>
void assign(MatrixRef<T, BaseType>& A, MatrixRef<T, BaseType> const& B)
{
   matrix_copy(B.as_derived()(), A.as_derived());
}

template <typename T, typename BaseType>
void add(MatrixRef<T, BaseType>& A, MatrixRef<T, BaseType> const& B)
{
   matrix_add(B.as_derived(), A.as_derived());
}

template <typename T, typename BaseType>
void subtract(MatrixRef<T, BaseType>& A, MatrixRef<T, BaseType> const& B)
{
   matrix_add_scaled(-number_traits<T>::identity(), B.as_derived(), A.as_derived());
}

// expression template for alpha * op(A)

template <typename T, typename BaseType>
class ScaledMatrix : public MatrixRef<T, BaseType, ScaledMatrix<T,BaseType>>
{
   public:
      using base_type = BaseType;

      ScaledMatrix(T const& Factor_, BaseType const& A_) : Factor(Factor_), A(A_) {}

      int rows() const { return A.rows(); }
      int cols() const { return A.cols(); }

      base_type const& base() const { return A; }
      T factor() const { return Factor; }

   private:
      T Factor;
      base_type const& A;
};

template <typename T, typename U, typename Derived, typename X>
ScaledMatrix<decltype(safe_convert<T>(std::declval<X>())), Derived>
operator*(X const& alpha, MatrixRef<T, U, Derived> const& M)
{
   return ScaledMatrix<T, Derived>(safe_convert<T>(alpha), M.as_derived());
}

template <typename T, typename BaseType>
void assign(MatrixRef<T, BaseType>& A, ScaledMatrix<T, BaseType> const& B)
{
   matrix_copy_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename BaseType>
void add(MatrixRef<T, BaseType>& A, ScaledMatrix<T, BaseType> const& B)
{
   matrix_add_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename BaseType>
void subtract(MatrixRef<T, BaseType>& A, ScaledMatrix<T, BaseType> const& B)
{
   matrix_add_scaled(-B.factor(), B.base(), A.as_derived());
}

// expression template for alpha * op(A) * op(B)

template <typename T, typename BaseType, typename U, typename V>
struct MatrixProduct : public MatrixRef<T, BaseType, MatrixProduct<T, BaseType, U, V>>
{
   MatrixProduct(MatrixRef<T, BaseType, U> const& A_, MatrixRef<T, BaseType, V> const& B_)
      : Factor(number_traits<T>::identity()), A(A_.as_derived()), B(B_.as_derived()) {}

   MatrixProduct(T const& Factor_, MatrixRef<T, BaseType, U> const& A_, MatrixRef<T, BaseType, V> const& B_)
      : Factor(Factor_), A(A_.as_derived()), B(B_.as_derived()) {}

   int rows() const { return A.rows(); }
   int cols() const { return B.cols(); }
   T factor() const { return Factor; }

   T Factor;
   U const& A;
   V const& B;
};

template <typename T, typename BaseType, typename U, typename V>
MatrixProduct<T, BaseType, U, V>
operator*(MatrixRef<T, BaseType, U> const& A, MatrixRef<T, BaseType, V> const& B)
{
   return MatrixProduct<T, BaseType, U, V>(A.as_derived(), B.as_derived());
}

template <typename T, typename BaseType, typename V>
MatrixProduct<T, BaseType, BaseType, V>
operator*(ScaledMatrix<T, BaseType> const& A, MatrixRef<T, BaseType, V> const& B)
{
   return MatrixProduct<T, BaseType, BaseType, V>(A.factor(), A.base(), B.as_derived());
}

// matrix product with a scalar, alpha * A * B

template <typename T, typename BaseType, typename Derived, typename U, typename V>
void assign(MatrixRef<T, BaseType, Derived>& C, MatrixProduct<T, BaseType, U, V> const& a)
{
   gemm(a.Factor, a.A, number_traits<T>::zero(), a.B, C.as_derived());
}

template <typename T, typename BaseType, typename Derived, typename U, typename V>
void add(MatrixRef<T, BaseType, Derived>& C, MatrixProduct<T, BaseType, U, V> const& a)
{
   gemm(a.Factor, a.A, number_traits<T>::identity(), a.B, C.as_derived());
}

template <typename T, typename BaseType, typename Derived, typename U, typename V>
void subtract(MatrixRef<T, BaseType, Derived>& C, MatrixProduct<T, BaseType, U, V> const& a)
{
   gemm(-a.Factor, a.A, number_traits<T>::identity(), a.B, C.as_derived());
}

} // namespace blas

#endif
