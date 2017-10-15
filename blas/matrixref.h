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
// DerivedType is the ultimate derived type, using the CRTP.
// Tag is a label that identifies storage type; matrices with the same
// storage can be freely mixed in expressions.  Typical use is
// for CPU versus GPU storage.
//
// The tag is also used for blas_traits<tag>.
//

template <typename ValueType, typename DerivedType, typename Tag>
class MatrixRef
{
   public:
      using value_type     = ValueType;
      using derived_type   = DerivedType;
      using tag_type       = Tag;

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

template <typename ValueType, typename DerivedType, typename Tag>
class BlasMatrix : public MatrixRef<ValueType, DerivedType, Tag>
{
   public:
      using value_type         = ValueType;
      using derived_type       = DerivedType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

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

template <typename T, typename BaseType, typename Tag>
class BlasMatrixTrans : public BlasMatrix<T, BlasMatrixTrans<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

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

template <typename T, typename BaseType, typename Tag>
class BlasMatrixHerm : public BlasMatrix<T, BlasMatrixHerm<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

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

template <typename T, typename BaseType, typename Tag>
class BlasMatrixConj : public BlasMatrix<T, BlasMatrixConj<T, BaseType, Tag>, Tag>
{
   public:
      using value_type         = T;
      using base_type          = BaseType;
      using tag_type           = Tag;
      using storage_type       = typename blas_traits<tag_type>::template storage_type<value_type>;
      using const_storage_type = typename blas_traits<tag_type>::template const_storage_type<value_type>;

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

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixTrans<T, BaseType, Tag>
trans(BlasMatrix<T, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<T, BaseType, Tag>(x.as_derived());
}

template <typename T, typename BaseType, typename Tag>
BaseType const&
trans(BlasMatrixTrans<T, BaseType, Tag> const& x)
{
   return x.base();
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixConj<T, BaseType, Tag>
trans(BlasMatrixHerm<T, BaseType, Tag> const& x)
{
   return BlasMatrixConj<T, BaseType, Tag>(x.base());
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixHerm<T, BaseType, Tag>
trans(BlasMatrixConj<T, BaseType, Tag> const& x)
{
   return BlasMatrixHerm<T, BaseType, Tag>(x.base());
}

// herm

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixHerm<T, BaseType, Tag>
herm(BlasMatrix<T, BaseType, Tag> const& x)
{
   return BlasMatrixHerm<T, BaseType, Tag>(x.as_derived());
}


template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixHerm<T, BaseType, Tag>
conj(BlasMatrixTrans<T, BaseType, Tag> const& x)
{
   return BlasMatrixHerm<T, BaseType, Tag>(x.base());
}

template <typename T, typename BaseType, typename Tag>
inline
BaseType const&
herm(BlasMatrixHerm<T, BaseType, Tag> const& x)
{
   return x.base();
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixTrans<T, BaseType, Tag>
herm(BlasMatrixConj<T, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<T, BaseType, Tag>(x.base());
}

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixConj<T, BaseType, Tag>
herm(BlasMatrixTrans<T, BaseType, Tag> const& x)
{
   return BlasMatrixConj<T, BaseType, Tag>(x.base());
}

// specialization for real types, herm() reduces to trans()
template <typename BaseType, typename Tag>
inline
BlasMatrixTrans<float, BaseType, Tag>
herm(BlasMatrix<float, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<float, BaseType, Tag>(x.as_derived());
}

template <typename BaseType, typename Tag>
inline
BlasMatrixTrans<double, BaseType, Tag>
herm(BlasMatrix<double, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<double, BaseType, Tag>(x.as_derived());
}

#if defined(HAVE_FLOAT128)
template <typename BaseType, typename Tag>
inline
BlasMatrixTrans<float128, BaseType, Tag>
herm(BlasMatrix<float128, BaseType, Tag> const& x)
{
   return BlasMatrixTrans<float128, BaseType, Tag>(x.as_derived());
}
#endif

// conj

template <typename T, typename BaseType, typename Tag>
inline
BlasMatrixConj<T, BaseType, Tag>
conj(BlasMatrix<T, BaseType, Tag> const& x)
{
   return BlasMatrixConj<T, BaseType, Tag>(x.as_derived());
}

// specialization for real types, conj() reduces to a no-op
template <typename BaseType, typename Tag>
inline
BaseType const&
conj(BlasMatrix<float, BaseType, Tag> const& x)
{
   return x.as_derived();
}

template <typename BaseType, typename Tag>
inline
BaseType const&
conj(BlasMatrix<double, BaseType, Tag> const& x)
{
   return x.as_derived();
}

#if defined(HAVE_FLOAT128)
template <typename BaseType, typename Tag>
inline
BaseType const&
conj(BlasMatrix<float128, BaseType, Tag> const& x)
{
   return x.as_derived();
}
#endif

// assignment

template <typename T, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_copy(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add(B.as_derived(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, U, Tag>& A, MatrixRef<T, V, Tag> const& B)
{
   matrix_add_scaled(-number_traits<T>::identity(), B.as_derived(), A.as_derived());
}

// expression template for alpha * op(A)

template <typename T, typename BaseType, typename Tag>
class ScaledMatrix : public MatrixRef<T, ScaledMatrix<T, BaseType, Tag>, Tag>
{
   public:
      using base_type = BaseType;

      ScaledMatrix(T const& Factor_, base_type const& A_) : Factor(Factor_), A(A_) {}

      int rows() const { return A.rows(); }
      int cols() const { return A.cols(); }

      base_type const& base() const { return A; }
      T factor() const { return Factor; }

   private:
      T Factor;
      base_type const& A;
};

template <typename T, typename BaseType, typename Tag, typename X>
inline
ScaledMatrix<decltype(safe_convert<T>(std::declval<X>())), BaseType, Tag>
operator*(X const& alpha, MatrixRef<T, BaseType, Tag> const& M)
{
   return ScaledMatrix<T, BaseType, Tag>(safe_convert<T>(alpha), M.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, U, Tag>& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_copy_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, U, Tag>& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(B.factor(), B.base(), A.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, U, Tag>& A, ScaledMatrix<T, V, Tag> const& B)
{
   matrix_add_scaled(-B.factor(), B.base(), A.as_derived());
}

// expression template for alpha * op(A) * op(B)

template <typename T, typename U, typename V, typename Tag>
struct MatrixProduct : public MatrixRef<T, MatrixProduct<T, U, V, Tag>, Tag>
{
   MatrixProduct(MatrixRef<T, U, Tag> const& A_, MatrixRef<T, V, Tag> const& B_)
      : Factor(number_traits<T>::identity()), A(A_.as_derived()), B(B_.as_derived()) {}

   MatrixProduct(T const& Factor_, MatrixRef<T, U, Tag> const& A_, MatrixRef<T, V, Tag> const& B_)
      : Factor(Factor_), A(A_.as_derived()), B(B_.as_derived()) {}

   int rows() const { return A.rows(); }
   int cols() const { return B.cols(); }
   T factor() const { return Factor; }

   T Factor;
   U const& A;
   V const& B;
};

template <typename T, typename U, typename V, typename Tag>
inline
MatrixProduct<T, U, V, Tag>
operator*(MatrixRef<T, U, Tag> const& A, MatrixRef<T, V, Tag> const& B)
{
   return MatrixProduct<T, U, V, Tag>(A.as_derived(), B.as_derived());
}

template <typename T, typename U, typename V, typename Tag>
inline
MatrixProduct<T, U, V, Tag>
operator*(ScaledMatrix<T, U, Tag> const& A, MatrixRef<T, V, Tag> const& B)
{
   return MatrixProduct<T, U, V, Tag>(A.factor(), A.base(), B.as_derived());
}

// matrix product with a scalar, alpha * A * B

template <typename T, typename Derived, typename U, typename V, typename Tag>
inline
void assign(MatrixRef<T, Derived, Tag>& C, MatrixProduct<T, U, V, Tag> const& a)
{
   gemm(a.Factor, a.A, number_traits<T>::zero(), a.B, C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
inline
void add(MatrixRef<T, Derived, Tag>& C, MatrixProduct<T, U, V, Tag> const& a)
{
   gemm(a.Factor, a.A, number_traits<T>::identity(), a.B, C.as_derived());
}

template <typename T, typename Derived, typename U, typename V, typename Tag>
inline
void subtract(MatrixRef<T, Derived, Tag>& C, MatrixProduct<T, U, V, Tag> const& a)
{
   gemm(-a.Factor, a.A, number_traits<T>::zero(), a.B, C.as_derived());
}

// miscellaneous

template <typename T, typename U, typename Tag>
void
inline
trace(BlasMatrix<T, U, Tag> const& x, typename blas_traits<Tag>::template async_ref<T>& r)
{
   vector_sum(x.as_derived().diagonal(), r);
}

template <typename T, typename U, typename Tag>
inline
T
trace(BlasMatrix<T, U, Tag> const& x)
{
   return vector_sum(x.as_derived().diagonal());
}

} // namespace blas

#endif
